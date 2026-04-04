#!/usr/bin/env python3
"""
Downloader CAMS global reanalysis (EAC4) per area Europa.

Comportamento chiave:
- Scarica e salva subito un file per ogni (profilo, variabile, mese).
- Supporta resume: i file gia' presenti vengono saltati.
- Supporta CAMS_VARIABLES=all con autodiscovery dal catalogo ADS.
- Supporta profili: single-level, model-level, pressure-level, oppure all.

Uso:
    python cams_eac4_download_europe.py
    python cams_eac4_download_europe.py --env .env
"""

from __future__ import annotations

import argparse
import calendar
import csv
import json
import logging
import os
import sys
import time
import urllib.request
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import cdsapi


DEFAULT_DATASET = "cams-global-reanalysis-eac4"
DEFAULT_API_URL = "https://ads.atmosphere.copernicus.eu/api"
DEFAULT_CATALOGUE_BASE_URL = "https://ads.atmosphere.copernicus.eu"
DEFAULT_TIME_STEPS = ["00:00", "03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00"]

STATUS_DOWNLOADED = "downloaded"
STATUS_SKIPPED = "skipped"
STATUS_ERROR = "error"
STATUS_DRY_RUN = "dry_run"

PROFILE_SINGLE = "single"
PROFILE_MODEL = "model"
PROFILE_PRESSURE = "pressure"


@dataclass
class ProfileSpec:
    name: str
    level_key: Optional[str]
    available_variables: Set[str]
    available_times: List[str]
    available_levels: List[str]
    available_date_range: Optional[Tuple[date, date]]


@dataclass
class ChunkResult:
    timestamp: str
    profile: str
    variable: str
    year: int
    month: int
    status: str
    throttle_detected: bool
    file_path: Path
    file_size_bytes: int
    note: str


@dataclass
class Config:
    dataset: str
    start_year: int
    start_month: int
    end_year: int
    end_month: int
    area: List[float]
    output_dir: Path
    file_format: str
    variables_raw: str
    level_type: str
    model_levels: List[str]
    pressure_levels: List[str]
    time_steps: List[str]
    max_retries: int
    retry_wait_seconds: int
    request_pause_seconds: float
    skip_existing: bool
    dry_run: bool
    cams_api_url: Optional[str]
    cams_api_key: Optional[str]
    catalogue_base_url: str
    log_file: Path
    log_level: str
    verify_download_magic: bool
    report_file: Path
    throttle_retry_wait_seconds: int
    throttle_max_wait_seconds: int
    adaptive_pause_enabled: bool
    adaptive_pause_increment_seconds: float
    adaptive_pause_decrement_seconds: float
    adaptive_pause_max_seconds: float


def parse_bool(value: Optional[str], default: bool) -> bool:
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "y", "on"}


def load_env_file(env_path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    if not env_path.exists():
        raise FileNotFoundError(f"File .env non trovato: {env_path}")

    for raw_line in env_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        if line.startswith("export "):
            line = line[7:].strip()

        if "=" not in line:
            continue

        key, raw_value = line.split("=", 1)
        key = key.strip()
        value = raw_value.strip().strip("\"'")
        values[key] = value

    return values


def get_env(values: Dict[str, str], key: str, default: Optional[str] = None) -> Optional[str]:
    env_val = os.getenv(key)
    if env_val is not None:
        return env_val
    return values.get(key, default)


def parse_list(raw: Optional[str]) -> List[str]:
    if not raw:
        return []
    parts = []
    for item in raw.replace("\n", ",").split(","):
        token = item.strip()
        if token:
            parts.append(token)
    return parts


def parse_area(raw: Optional[str]) -> List[float]:
    if not raw:
        raise ValueError("Manca CAMS_AREA nel file .env")

    parts = [p.strip() for p in raw.split(",") if p.strip()]
    if len(parts) != 4:
        raise ValueError("CAMS_AREA deve avere 4 valori: nord,ovest,sud,est (es: 72,-25,34,45)")

    area = [float(p) for p in parts]
    north, west, south, east = area
    if north <= south:
        raise ValueError("CAMS_AREA non valida: north deve essere > south")
    if east <= west:
        raise ValueError("CAMS_AREA non valida: east deve essere > west")
    return area


def parse_levels(raw: Optional[str]) -> List[str]:
    values = parse_list(raw)
    return sorted(set(values), key=lambda x: int(x)) if values else []


def parse_time_steps(raw: Optional[str]) -> List[str]:
    values = parse_list(raw)
    if not values:
        return DEFAULT_TIME_STEPS.copy()

    normalized: List[str] = []
    for item in values:
        token = item
        if len(token) == 2 and token.isdigit():
            token = f"{token}:00"
        if len(token) != 5 or token[2] != ":":
            raise ValueError(f"CAMS_TIME_STEPS non valido: {item}")
        hh = token[:2]
        mm = token[3:]
        if not (hh.isdigit() and mm.isdigit()):
            raise ValueError(f"CAMS_TIME_STEPS non valido: {item}")
        h = int(hh)
        m = int(mm)
        if h < 0 or h > 23 or m < 0 or m > 59:
            raise ValueError(f"CAMS_TIME_STEPS non valido: {item}")
        normalized.append(f"{h:02d}:{m:02d}")

    return sorted(set(normalized))


def parse_date_range(raw: str) -> Tuple[date, date]:
    # ADS constraints usa formato "YYYY-MM-DD/YYYY-MM-DD"
    start_raw, end_raw = raw.split("/", 1)
    start = date.fromisoformat(start_raw)
    end = date.fromisoformat(end_raw)
    return start, end


def parse_config(env_values: Dict[str, str]) -> Config:
    start_year = int(get_env(env_values, "CAMS_START_YEAR", "2003") or "2003")
    start_month = int(get_env(env_values, "CAMS_START_MONTH", "1") or "1")
    end_year = int(get_env(env_values, "CAMS_END_YEAR", "2025") or "2025")
    end_month = int(get_env(env_values, "CAMS_END_MONTH", "12") or "12")

    if not (1 <= start_month <= 12 and 1 <= end_month <= 12):
        raise ValueError("I mesi devono essere compresi tra 1 e 12")

    start = date(start_year, start_month, 1)
    end = date(end_year, end_month, 1)
    if start > end:
        raise ValueError("Intervallo non valido: data di inizio successiva alla fine")

    output_dir_raw = get_env(env_values, "CAMS_OUTPUT_DIR")
    if not output_dir_raw:
        raise ValueError("Manca CAMS_OUTPUT_DIR nel file .env")

    output_dir = Path(output_dir_raw).expanduser()

    file_format = (get_env(env_values, "CAMS_FILE_FORMAT", "grib") or "grib").strip().lower()
    if file_format == "netcdf":
        file_format = "netcdf_zip"
    if file_format not in {"grib", "netcdf_zip"}:
        raise ValueError("CAMS_FILE_FORMAT deve essere 'grib', 'netcdf' oppure 'netcdf_zip'")

    level_type = (get_env(env_values, "CAMS_LEVEL_TYPE", "single") or "single").strip().lower()
    if level_type not in {PROFILE_SINGLE, PROFILE_MODEL, PROFILE_PRESSURE, "all"}:
        raise ValueError("CAMS_LEVEL_TYPE deve essere: single, model, pressure, all")

    log_file_raw = get_env(env_values, "CAMS_LOG_FILE", "") or ""
    log_file = Path(log_file_raw).expanduser() if log_file_raw.strip() else output_dir / "cams_eac4_downloader.log"

    report_file_raw = get_env(env_values, "CAMS_REPORT_FILE", "") or ""
    report_file = Path(report_file_raw).expanduser() if report_file_raw.strip() else output_dir / "cams_eac4_report.csv"

    # Priorita': CAMSAPI_*, poi CDSAPI_*
    cams_api_url = get_env(env_values, "CAMSAPI_URL") or get_env(env_values, "CDSAPI_URL")
    cams_api_key = get_env(env_values, "CAMSAPI_KEY") or get_env(env_values, "CDSAPI_KEY")

    return Config(
        dataset=get_env(env_values, "CAMS_DATASET", DEFAULT_DATASET) or DEFAULT_DATASET,
        start_year=start_year,
        start_month=start_month,
        end_year=end_year,
        end_month=end_month,
        area=parse_area(get_env(env_values, "CAMS_AREA", "72,-25,34,45")),
        output_dir=output_dir,
        file_format=file_format,
        variables_raw=(get_env(env_values, "CAMS_VARIABLES", "all") or "all").strip(),
        level_type=level_type,
        model_levels=parse_levels(get_env(env_values, "CAMS_MODEL_LEVELS", "")),
        pressure_levels=parse_levels(get_env(env_values, "CAMS_PRESSURE_LEVELS", "")),
        time_steps=parse_time_steps(get_env(env_values, "CAMS_TIME_STEPS", "")),
        max_retries=int(get_env(env_values, "CAMS_MAX_RETRIES", "3") or "3"),
        retry_wait_seconds=int(get_env(env_values, "CAMS_RETRY_WAIT_SECONDS", "25") or "25"),
        request_pause_seconds=float(get_env(env_values, "CAMS_REQUEST_PAUSE_SECONDS", "2") or "2"),
        skip_existing=parse_bool(get_env(env_values, "CAMS_SKIP_EXISTING", "true"), True),
        dry_run=parse_bool(get_env(env_values, "CAMS_DRY_RUN", "false"), False),
        cams_api_url=cams_api_url,
        cams_api_key=cams_api_key,
        catalogue_base_url=get_env(env_values, "CAMS_CATALOGUE_BASE_URL", DEFAULT_CATALOGUE_BASE_URL)
        or DEFAULT_CATALOGUE_BASE_URL,
        log_file=log_file,
        log_level=(get_env(env_values, "CAMS_LOG_LEVEL", "INFO") or "INFO").upper(),
        verify_download_magic=parse_bool(get_env(env_values, "CAMS_VERIFY_DOWNLOAD_MAGIC", "true"), True),
        report_file=report_file,
        throttle_retry_wait_seconds=int(get_env(env_values, "CAMS_THROTTLE_RETRY_WAIT_SECONDS", "120") or "120"),
        throttle_max_wait_seconds=int(get_env(env_values, "CAMS_THROTTLE_MAX_WAIT_SECONDS", "900") or "900"),
        adaptive_pause_enabled=parse_bool(get_env(env_values, "CAMS_ADAPTIVE_PAUSE_ENABLED", "true"), True),
        adaptive_pause_increment_seconds=float(get_env(env_values, "CAMS_ADAPTIVE_PAUSE_INCREMENT_SECONDS", "2") or "2"),
        adaptive_pause_decrement_seconds=float(get_env(env_values, "CAMS_ADAPTIVE_PAUSE_DECREMENT_SECONDS", "0.5") or "0.5"),
        adaptive_pause_max_seconds=float(get_env(env_values, "CAMS_ADAPTIVE_PAUSE_MAX_SECONDS", "30") or "30"),
    )


def month_iter(start_year: int, start_month: int, end_year: int, end_month: int) -> Iterable[Tuple[int, int]]:
    year = start_year
    month = start_month
    while (year < end_year) or (year == end_year and month <= end_month):
        yield year, month
        month += 1
        if month > 12:
            month = 1
            year += 1


def month_date_range(year: int, month: int) -> str:
    last_day = calendar.monthrange(year, month)[1]
    return f"{year:04d}-{month:02d}-01/{year:04d}-{month:02d}-{last_day:02d}"


def build_target_path(output_dir: Path, profile: str, variable: str, year: int, month: int, file_format: str) -> Path:
    ext = "grib" if file_format == "grib" else "zip"
    return output_dir / profile / variable / f"cams_eac4_{profile}_{variable}_{year}{month:02d}.{ext}"


def _http_get_json(url: str) -> object:
    req = urllib.request.Request(url=url, headers={"User-Agent": "cams-eac4-downloader/1.0"})
    with urllib.request.urlopen(req, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_profile_specs(config: Config) -> Dict[str, ProfileSpec]:
    collection_url = f"{config.catalogue_base_url.rstrip('/')}/api/catalogue/v1/collections/{config.dataset}"
    collection = _http_get_json(collection_url)
    if not isinstance(collection, dict):
        raise RuntimeError("Metadata collection non valida dal catalogo ADS")

    links = collection.get("links") or []
    constraints_url: Optional[str] = None
    for link in links:
        if isinstance(link, dict) and link.get("rel") == "constraints":
            constraints_url = link.get("href")
            break

    if not constraints_url:
        raise RuntimeError("Impossibile trovare link constraints nel catalogo ADS")

    constraints = _http_get_json(constraints_url)
    if not isinstance(constraints, list):
        raise RuntimeError("Formato constraints non valido")

    specs: Dict[str, ProfileSpec] = {
        PROFILE_SINGLE: ProfileSpec(
            name=PROFILE_SINGLE,
            level_key=None,
            available_variables=set(),
            available_times=[],
            available_levels=[],
            available_date_range=None,
        ),
        PROFILE_MODEL: ProfileSpec(
            name=PROFILE_MODEL,
            level_key="model_level",
            available_variables=set(),
            available_times=[],
            available_levels=[],
            available_date_range=None,
        ),
        PROFILE_PRESSURE: ProfileSpec(
            name=PROFILE_PRESSURE,
            level_key="pressure_level",
            available_variables=set(),
            available_times=[],
            available_levels=[],
            available_date_range=None,
        ),
    }

    for item in constraints:
        if not isinstance(item, dict):
            continue

        if "pressure_level" in item:
            profile_name = PROFILE_PRESSURE
            level_key = "pressure_level"
        elif "model_level" in item:
            profile_name = PROFILE_MODEL
            level_key = "model_level"
        else:
            profile_name = PROFILE_SINGLE
            level_key = None

        spec = specs[profile_name]

        raw_vars = item.get("variable") or []
        if isinstance(raw_vars, list):
            spec.available_variables.update(v for v in raw_vars if isinstance(v, str))

        raw_times = item.get("time") or []
        if isinstance(raw_times, list):
            for tkn in raw_times:
                if isinstance(tkn, str) and tkn not in spec.available_times:
                    spec.available_times.append(tkn)

        if level_key:
            raw_levels = item.get(level_key) or []
            if isinstance(raw_levels, list):
                for lvl in raw_levels:
                    if isinstance(lvl, str) and lvl not in spec.available_levels:
                        spec.available_levels.append(lvl)

        raw_dates = item.get("date") or []
        if isinstance(raw_dates, list) and raw_dates:
            try:
                start, end = parse_date_range(str(raw_dates[0]))
                if spec.available_date_range is None:
                    spec.available_date_range = (start, end)
                else:
                    old_start, old_end = spec.available_date_range
                    spec.available_date_range = (min(old_start, start), max(old_end, end))
            except Exception:
                pass

    for spec in specs.values():
        spec.available_times.sort()
        spec.available_levels.sort(key=lambda x: int(x))

    return specs


def get_client(config: Config) -> cdsapi.Client:
    kwargs: Dict[str, object] = {"quiet": False, "progress": True}
    if config.cams_api_key:
        kwargs["url"] = config.cams_api_url or DEFAULT_API_URL
        kwargs["key"] = config.cams_api_key
    elif config.cams_api_url:
        kwargs["url"] = config.cams_api_url
    return cdsapi.Client(**kwargs)


def setup_logger(config: Config) -> logging.Logger:
    logger = logging.getLogger("cams_eac4_downloader")
    logger.handlers.clear()

    level = getattr(logging, config.log_level, logging.INFO)
    logger.setLevel(level)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(level)
    stream_handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))
    logger.addHandler(stream_handler)

    config.log_file.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(config.log_file, encoding="utf-8")
    file_handler.setLevel(level)
    file_handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))
    logger.addHandler(file_handler)

    return logger


def validate_target_path(target: Path, output_dir: Path, file_format: str) -> Tuple[bool, str]:
    expected_suffix = ".grib" if file_format == "grib" else ".zip"
    if target.suffix.lower() != expected_suffix:
        return False, f"Estensione non coerente con file_format: {target.suffix} != {expected_suffix}"

    output_abs = output_dir.resolve()
    target_abs = target.resolve()
    try:
        target_abs.relative_to(output_abs)
    except ValueError:
        return False, f"Target fuori da CAMS_OUTPUT_DIR: {target_abs}"
    return True, "ok"


def verify_downloaded_file(target: Path, file_format: str, verify_magic: bool) -> Tuple[bool, str]:
    if not target.exists():
        return False, "file non trovato dopo il download"

    size = target.stat().st_size
    if size <= 0:
        return False, "file vuoto"

    if not verify_magic:
        return True, f"ok, size={size}"

    with target.open("rb") as fh:
        magic = fh.read(8)

    if file_format == "grib":
        if not magic.startswith(b"GRIB"):
            return False, f"header GRIB non valido: {magic!r}"
    else:
        # netcdf_zip viene consegnato come archivio zip
        if not magic.startswith(b"PK"):
            return False, f"header ZIP non valido per netcdf_zip: {magic!r}"

    return True, f"ok, size={size}"


def is_throttle_error(exc: Exception) -> bool:
    msg = str(exc).lower()
    throttle_patterns = ["429", "too many requests", "rate limit", "rate-limit", "throttl"]
    return any(p in msg for p in throttle_patterns)


def validate_date_against_profile(
    profile: ProfileSpec,
    year: int,
    month: int,
) -> Tuple[bool, str]:
    if profile.available_date_range is None:
        return True, "ok"

    first = date(year, month, 1)
    last = date(year, month, calendar.monthrange(year, month)[1])
    start, end = profile.available_date_range
    if first < start or last > end:
        return False, f"periodo fuori copertura profilo ({start.isoformat()} -> {end.isoformat()})"
    return True, "ok"


def select_profiles(config: Config, specs: Dict[str, ProfileSpec]) -> List[ProfileSpec]:
    if config.level_type == "all":
        selected = [specs[PROFILE_SINGLE], specs[PROFILE_MODEL], specs[PROFILE_PRESSURE]]
    else:
        selected = [specs[config.level_type]]

    return [p for p in selected if p.available_variables]


def select_variables_for_profile(config: Config, profile: ProfileSpec) -> List[str]:
    if config.variables_raw.lower() == "all":
        return sorted(profile.available_variables)

    requested = parse_list(config.variables_raw)
    if not requested:
        raise ValueError("CAMS_VARIABLES non valido: inserire almeno una variabile o 'all'")

    missing = [v for v in requested if v not in profile.available_variables]
    if missing:
        raise ValueError(
            f"Variabili non disponibili per profilo '{profile.name}': {', '.join(sorted(missing))}"
        )

    return sorted(set(requested))


def select_levels_for_profile(config: Config, profile: ProfileSpec) -> List[str]:
    if profile.level_key is None:
        return []

    if profile.level_key == "model_level":
        requested = config.model_levels
    else:
        requested = config.pressure_levels

    if not requested:
        return profile.available_levels.copy()

    missing = [lvl for lvl in requested if lvl not in profile.available_levels]
    if missing:
        raise ValueError(
            f"Livelli non disponibili per {profile.level_key}: {', '.join(sorted(missing, key=lambda x: int(x)))}"
        )

    return requested


def select_times_for_profile(config: Config, profile: ProfileSpec) -> List[str]:
    if not profile.available_times:
        return config.time_steps

    missing = [t for t in config.time_steps if t not in profile.available_times]
    if missing:
        raise ValueError(
            f"Time step non disponibili nel profilo '{profile.name}': {', '.join(sorted(missing))}"
        )

    return config.time_steps


def build_request(
    config: Config,
    profile: ProfileSpec,
    variable: str,
    year: int,
    month: int,
    profile_levels: List[str],
    profile_times: List[str],
) -> Dict[str, object]:
    request: Dict[str, object] = {
        "variable": [variable],
        "date": month_date_range(year, month),
        "time": profile_times,
        "area": config.area,
        "data_format": config.file_format,
    }
    if profile.level_key:
        request[profile.level_key] = profile_levels
    return request


def download_single_chunk(
    client: cdsapi.Client,
    config: Config,
    logger: logging.Logger,
    profile: ProfileSpec,
    variable: str,
    year: int,
    month: int,
    profile_levels: List[str],
    profile_times: List[str],
) -> ChunkResult:
    ts = datetime.now().isoformat(timespec="seconds")
    target = build_target_path(config.output_dir, profile.name, variable, year, month, config.file_format)
    target.parent.mkdir(parents=True, exist_ok=True)

    def _result(status: str, note: str, throttle_detected: bool = False) -> ChunkResult:
        size = target.stat().st_size if target.exists() else 0
        return ChunkResult(
            timestamp=ts,
            profile=profile.name,
            variable=variable,
            year=year,
            month=month,
            status=status,
            throttle_detected=throttle_detected,
            file_path=target,
            file_size_bytes=size,
            note=note,
        )

    date_ok, date_msg = validate_date_against_profile(profile, year, month)
    if not date_ok:
        logger.error(
            "[DATE_INVALID] profile=%s variable=%s year=%s month=%s msg=%s",
            profile.name,
            variable,
            year,
            month,
            date_msg,
        )
        return _result(STATUS_ERROR, f"DATE_INVALID: {date_msg}")

    target_ok, target_msg = validate_target_path(target, config.output_dir, config.file_format)
    if not target_ok:
        logger.error(
            "[TARGET_INVALID] profile=%s variable=%s year=%s month=%s target=%s msg=%s",
            profile.name,
            variable,
            year,
            month,
            target,
            target_msg,
        )
        return _result(STATUS_ERROR, f"TARGET_INVALID: {target_msg}")

    request = build_request(config, profile, variable, year, month, profile_levels, profile_times)

    logger.info(
        "[REQUEST_READY] dataset=%s profile=%s variable=%s year=%s month=%02d format=%s target=%s",
        config.dataset,
        profile.name,
        variable,
        year,
        month,
        config.file_format,
        target,
    )

    if config.skip_existing and target.exists() and target.stat().st_size > 0:
        ok, msg = verify_downloaded_file(target, config.file_format, config.verify_download_magic)
        if ok:
            logger.info(
                "[SKIP_EXISTING_OK] profile=%s variable=%s year=%s month=%02d target=%s %s",
                profile.name,
                variable,
                year,
                month,
                target,
                msg,
            )
            return _result(STATUS_SKIPPED, msg)
        logger.warning(
            "[SKIP_EXISTING_BAD] profile=%s variable=%s year=%s month=%02d target=%s msg=%s",
            profile.name,
            variable,
            year,
            month,
            target,
            msg,
        )
        target.unlink(missing_ok=True)

    if config.dry_run:
        logger.info(
            "[DRY_RUN] profile=%s variable=%s year=%s month=%02d target=%s",
            profile.name,
            variable,
            year,
            month,
            target,
        )
        return _result(STATUS_DRY_RUN, "dry_run")

    attempt = 0
    last_error = ""
    throttle_detected = False
    while attempt < config.max_retries:
        attempt += 1
        try:
            logger.info(
                "[DOWNLOAD_START] profile=%s variable=%s year=%s month=%02d attempt=%s/%s",
                profile.name,
                variable,
                year,
                month,
                attempt,
                config.max_retries,
            )
            client.retrieve(config.dataset, request, str(target))
            ok, msg = verify_downloaded_file(target, config.file_format, config.verify_download_magic)
            if not ok:
                logger.error(
                    "[VERIFY_FAILED] profile=%s variable=%s year=%s month=%02d target=%s msg=%s",
                    profile.name,
                    variable,
                    year,
                    month,
                    target,
                    msg,
                )
                target.unlink(missing_ok=True)
                raise RuntimeError(msg)

            logger.info(
                "[DOWNLOAD_DONE] profile=%s variable=%s year=%s month=%02d target=%s %s",
                profile.name,
                variable,
                year,
                month,
                target,
                msg,
            )
            return _result(STATUS_DOWNLOADED, msg, throttle_detected=throttle_detected)
        except Exception as exc:
            last_error = str(exc)
            current_is_throttle = is_throttle_error(exc)
            throttle_detected = throttle_detected or current_is_throttle

            logger.error(
                "[DOWNLOAD_ERROR] profile=%s variable=%s year=%s month=%02d attempt=%s/%s target=%s err=%s",
                profile.name,
                variable,
                year,
                month,
                attempt,
                config.max_retries,
                target,
                exc,
            )
            if current_is_throttle:
                logger.warning(
                    "[THROTTLE_DETECTED] profile=%s variable=%s year=%s month=%02d attempt=%s/%s",
                    profile.name,
                    variable,
                    year,
                    month,
                    attempt,
                    config.max_retries,
                )

            if target.exists() and target.stat().st_size == 0:
                target.unlink(missing_ok=True)

            if attempt < config.max_retries:
                if current_is_throttle:
                    throttle_wait = min(
                        config.throttle_max_wait_seconds,
                        max(config.throttle_retry_wait_seconds, config.retry_wait_seconds) * (2 ** (attempt - 1)),
                    )
                    logger.info(
                        "[THROTTLE_WAIT] profile=%s variable=%s year=%s month=%02d wait_seconds=%s",
                        profile.name,
                        variable,
                        year,
                        month,
                        throttle_wait,
                    )
                    time.sleep(throttle_wait)
                else:
                    logger.info(
                        "[RETRY_WAIT] profile=%s variable=%s year=%s month=%02d wait_seconds=%s",
                        profile.name,
                        variable,
                        year,
                        month,
                        config.retry_wait_seconds,
                    )
                    time.sleep(config.retry_wait_seconds)

    return _result(STATUS_ERROR, f"max retries esauriti: {last_error}", throttle_detected=throttle_detected)


def _write_report(config: Config, results: List[ChunkResult], logger: logging.Logger) -> None:
    report_path = config.report_file
    report_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not report_path.exists()

    with report_path.open("a", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        if write_header:
            writer.writerow(
                [
                    "timestamp",
                    "profile",
                    "variable",
                    "year",
                    "month",
                    "status",
                    "throttle_detected",
                    "file_path",
                    "file_size_bytes",
                    "note",
                ]
            )
        for r in results:
            writer.writerow(
                [
                    r.timestamp,
                    r.profile,
                    r.variable,
                    r.year,
                    f"{r.month:02d}",
                    r.status,
                    str(r.throttle_detected).lower(),
                    str(r.file_path),
                    r.file_size_bytes,
                    r.note,
                ]
            )

    logger.info("[REPORT] Scritto: %s (%s righe)", report_path, len(results))


def run(config: Config) -> int:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(config)

    logger.info("=" * 80)
    logger.info("CAMS EAC4 downloader")
    logger.info("Dataset      : %s", config.dataset)
    logger.info("Periodo      : %s-%02d -> %s-%02d", config.start_year, config.start_month, config.end_year, config.end_month)
    logger.info("Area         : N/W/S/E = %s", config.area)
    logger.info("Output       : %s", config.output_dir)
    logger.info("Formato      : %s", config.file_format)
    logger.info("Level type   : %s", config.level_type)
    logger.info("Variables    : %s", config.variables_raw)
    logger.info("Time steps   : %s", ",".join(config.time_steps))
    logger.info("Skip existing: %s", config.skip_existing)
    logger.info("Dry run      : %s", config.dry_run)
    logger.info("Log file     : %s", config.log_file)
    logger.info("Throttle wait: %ss (max %ss)", config.throttle_retry_wait_seconds, config.throttle_max_wait_seconds)
    logger.info("Adaptive pause: enabled=%s base=%ss max=%ss", config.adaptive_pause_enabled, config.request_pause_seconds, config.adaptive_pause_max_seconds)
    logger.info("=" * 80)

    profile_specs = fetch_profile_specs(config)
    selected_profiles = select_profiles(config, profile_specs)
    if not selected_profiles:
        logger.error("Nessun profilo disponibile per la selezione corrente")
        return 2

    profile_context: Dict[str, Dict[str, object]] = {}
    for profile in selected_profiles:
        variables = select_variables_for_profile(config, profile)
        levels = select_levels_for_profile(config, profile)
        times = select_times_for_profile(config, profile)
        profile_context[profile.name] = {
            "variables": variables,
            "levels": levels,
            "times": times,
        }
        logger.info(
            "[PROFILE] name=%s variables=%s levels=%s times=%s coverage=%s",
            profile.name,
            len(variables),
            len(levels) if profile.level_key else 0,
            len(times),
            (
                f"{profile.available_date_range[0].isoformat()}->{profile.available_date_range[1].isoformat()}"
                if profile.available_date_range
                else "unknown"
            ),
        )

    client = get_client(config)

    results: List[ChunkResult] = []
    current_pause_seconds = config.request_pause_seconds

    for profile in selected_profiles:
        ctx = profile_context[profile.name]
        profile_variables = ctx["variables"]
        profile_levels = ctx["levels"]
        profile_times = ctx["times"]

        for variable in profile_variables:
            for year, month in month_iter(
                config.start_year,
                config.start_month,
                config.end_year,
                config.end_month,
            ):
                result = download_single_chunk(
                    client=client,
                    config=config,
                    logger=logger,
                    profile=profile,
                    variable=variable,
                    year=year,
                    month=month,
                    profile_levels=profile_levels,
                    profile_times=profile_times,
                )
                results.append(result)

                if config.adaptive_pause_enabled:
                    if result.throttle_detected:
                        new_pause = min(
                            config.adaptive_pause_max_seconds,
                            current_pause_seconds + config.adaptive_pause_increment_seconds,
                        )
                        if new_pause > current_pause_seconds:
                            logger.info(
                                "[ADAPTIVE_PAUSE_UP] old=%ss new=%ss reason=throttle",
                                round(current_pause_seconds, 3),
                                round(new_pause, 3),
                            )
                            current_pause_seconds = new_pause
                    elif result.status in {STATUS_DOWNLOADED, STATUS_SKIPPED, STATUS_DRY_RUN} and current_pause_seconds > config.request_pause_seconds:
                        new_pause = max(
                            config.request_pause_seconds,
                            current_pause_seconds - config.adaptive_pause_decrement_seconds,
                        )
                        if new_pause < current_pause_seconds:
                            logger.info(
                                "[ADAPTIVE_PAUSE_DOWN] old=%ss new=%ss reason=no_throttle",
                                round(current_pause_seconds, 3),
                                round(new_pause, 3),
                            )
                            current_pause_seconds = new_pause

                if result.status != STATUS_SKIPPED and current_pause_seconds > 0:
                    logger.debug("[INTER_CHUNK_PAUSE] seconds=%s", round(current_pause_seconds, 3))
                    time.sleep(current_pause_seconds)

    n_downloaded = sum(1 for r in results if r.status == STATUS_DOWNLOADED)
    n_skipped = sum(1 for r in results if r.status == STATUS_SKIPPED)
    n_dry_run = sum(1 for r in results if r.status == STATUS_DRY_RUN)
    n_error = sum(1 for r in results if r.status == STATUS_ERROR)
    n_throttle = sum(1 for r in results if r.throttle_detected)
    total = len(results)

    logger.info("=" * 80)
    logger.info("Riepilogo finale:")
    logger.info("  Totale chunk pianificati : %s", total)
    logger.info("  Scaricati (nuovi)        : %s", n_downloaded)
    logger.info("  Saltati (gia' presenti)  : %s", n_skipped)
    logger.info("  Dry-run                  : %s", n_dry_run)
    logger.info("  Errori                   : %s", n_error)
    logger.info("  Chunk con throttling     : %s", n_throttle)
    logger.info("=" * 80)

    _write_report(config, results, logger)
    return 0 if n_error == 0 else 2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download CAMS global reanalysis EAC4 per Europa")
    parser.add_argument(
        "--env",
        default=".env",
        help="Percorso del file .env con i parametri (default: .env)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    env_path = Path(args.env).expanduser().resolve()

    try:
        env_values = load_env_file(env_path)
        config = parse_config(env_values)
    except Exception as exc:
        print(f"Errore configurazione: {exc}", file=sys.stderr)
        return 1

    try:
        return run(config)
    except Exception as exc:
        print(f"Errore runtime: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
