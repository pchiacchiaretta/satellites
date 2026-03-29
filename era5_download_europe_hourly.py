#!/usr/bin/env python3
"""
Downloader ERA5 hourly single levels (Copernicus CDS) per Europa.

Comportamento chiave:
- Scarica e salva subito un file per ogni (variabile, mese).
- Non aspetta la fine dell'intero job prima di scrivere su disco.
- Supporta resume: i file gia' presenti vengono saltati.

Prerequisiti:
- account CDS e API key attiva
- pacchetto Python: cdsapi

Uso:
    python era5_download_europe_hourly.py
    python era5_download_europe_hourly.py --env .env
"""

from __future__ import annotations

import argparse
import calendar
import csv
import logging
import os
import sys
import time
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import cdsapi


DEFAULT_VARIABLES = [
    "2m_temperature",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "convective_precipitation",
    "convective_rain_rate",
    "instantaneous_10m_wind_gust",
    "k_index",
    "precipitation_type",
    "total_cloud_cover",
    "total_column_ozone",
    "total_column_rain_water",
    "total_column_snow_water",
    "total_column_supercooled_liquid_water",
    "total_column_water",
    "total_column_water_vapour",
    "total_precipitation",
]

TIME_STEPS = [f"{hour:02d}:00" for hour in range(24)]
DEFAULT_CDS_API_URL = "https://cds.climate.copernicus.eu/api"

# Possibili stati di un singolo chunk di download
STATUS_DOWNLOADED = "downloaded"
STATUS_SKIPPED    = "skipped"
STATUS_ERROR      = "error"
STATUS_DRY_RUN    = "dry_run"


@dataclass
class ChunkResult:
    timestamp: str
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
    product_type: str
    start_year: int
    start_month: int
    end_year: int
    end_month: int
    area: List[float]
    output_dir: Path
    file_format: str
    variables: List[str]
    max_retries: int
    retry_wait_seconds: int
    request_pause_seconds: float
    skip_existing: bool
    dry_run: bool
    cds_url: str | None
    cds_key: str | None
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


def parse_bool(value: str, default: bool) -> bool:
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


def get_env(values: Dict[str, str], key: str, default: str | None = None) -> str | None:
    # Le variabili d'ambiente reali (es. Docker Compose environment:) hanno
    # priorità sul file .env, cosi' gli override del container funzionano.
    env_val = os.getenv(key)
    if env_val is not None:
        return env_val
    return values.get(key, default)


def parse_variables(raw: str | None) -> List[str]:
    if not raw:
        return DEFAULT_VARIABLES.copy()

    variables: List[str] = []
    for part in raw.replace("\n", ",").split(","):
        variable = part.strip()
        if variable:
            variables.append(variable)

    if not variables:
        raise ValueError("ERA5_VARIABLES e' vuoto: inserisci almeno una variabile.")
    return variables


def parse_area(raw: str | None) -> List[float]:
    if not raw:
        raise ValueError("Manca ERA5_AREA nel file .env")

    parts = [p.strip() for p in raw.split(",") if p.strip()]
    if len(parts) != 4:
        raise ValueError(
            "ERA5_AREA deve avere 4 valori: nord,ovest,sud,est (es: 72,-25,34,45)"
        )

    area = [float(p) for p in parts]
    north, west, south, east = area
    if north <= south:
        raise ValueError("ERA5_AREA non valida: north deve essere > south")
    if east <= west:
        raise ValueError("ERA5_AREA non valida: east deve essere > west")

    return area


def parse_config(env_values: Dict[str, str]) -> Config:
    start_year = int(get_env(env_values, "ERA5_START_YEAR", "2014"))
    start_month = int(get_env(env_values, "ERA5_START_MONTH", "1"))
    end_year = int(get_env(env_values, "ERA5_END_YEAR", "2026"))
    end_month = int(get_env(env_values, "ERA5_END_MONTH", "12"))

    if not (1 <= start_month <= 12 and 1 <= end_month <= 12):
        raise ValueError("I mesi devono essere compresi tra 1 e 12")

    start = date(start_year, start_month, 1)
    end = date(end_year, end_month, 1)
    if start > end:
        raise ValueError("Intervallo non valido: data di inizio successiva alla fine")

    output_dir_value = get_env(env_values, "ERA5_OUTPUT_DIR")
    if not output_dir_value:
        raise ValueError("Manca ERA5_OUTPUT_DIR nel file .env")

    file_format = (get_env(env_values, "ERA5_FILE_FORMAT", "grib") or "grib").lower().strip()
    if file_format not in {"grib", "netcdf"}:
        raise ValueError("ERA5_FILE_FORMAT deve essere 'grib' oppure 'netcdf'")

    output_dir_path = Path(output_dir_value).expanduser()
    log_file_raw = get_env(env_values, "ERA5_LOG_FILE", "") or ""
    log_file = Path(log_file_raw).expanduser() if log_file_raw.strip() else output_dir_path / "era5_downloader.log"

    report_file_raw = get_env(env_values, "ERA5_REPORT_FILE", "") or ""
    report_file = Path(report_file_raw).expanduser() if report_file_raw.strip() else output_dir_path / "era5_report.csv"

    return Config(
        dataset=get_env(env_values, "ERA5_DATASET", "reanalysis-era5-single-levels")
        or "reanalysis-era5-single-levels",
        product_type=get_env(env_values, "ERA5_PRODUCT_TYPE", "reanalysis")
        or "reanalysis",
        start_year=start_year,
        start_month=start_month,
        end_year=end_year,
        end_month=end_month,
        area=parse_area(get_env(env_values, "ERA5_AREA")),
        output_dir=output_dir_path,
        file_format=file_format,
        variables=parse_variables(get_env(env_values, "ERA5_VARIABLES")),
        max_retries=int(get_env(env_values, "ERA5_MAX_RETRIES", "3") or "3"),
        retry_wait_seconds=int(get_env(env_values, "ERA5_RETRY_WAIT_SECONDS", "25") or "25"),
        request_pause_seconds=float(get_env(env_values, "ERA5_REQUEST_PAUSE_SECONDS", "2") or "2"),
        skip_existing=parse_bool(get_env(env_values, "ERA5_SKIP_EXISTING", "true"), True),
        dry_run=parse_bool(get_env(env_values, "ERA5_DRY_RUN", "false"), False),
        cds_url=get_env(env_values, "CDSAPI_URL"),
        cds_key=get_env(env_values, "CDSAPI_KEY"),
        log_file=log_file,
        log_level=(get_env(env_values, "ERA5_LOG_LEVEL", "INFO") or "INFO").upper(),
        verify_download_magic=parse_bool(get_env(env_values, "ERA5_VERIFY_DOWNLOAD_MAGIC", "true"), True),
        report_file=report_file,
        throttle_retry_wait_seconds=int(get_env(env_values, "ERA5_THROTTLE_RETRY_WAIT_SECONDS", "120") or "120"),
        throttle_max_wait_seconds=int(get_env(env_values, "ERA5_THROTTLE_MAX_WAIT_SECONDS", "900") or "900"),
        adaptive_pause_enabled=parse_bool(get_env(env_values, "ERA5_ADAPTIVE_PAUSE_ENABLED", "true"), True),
        adaptive_pause_increment_seconds=float(get_env(env_values, "ERA5_ADAPTIVE_PAUSE_INCREMENT_SECONDS", "2") or "2"),
        adaptive_pause_decrement_seconds=float(get_env(env_values, "ERA5_ADAPTIVE_PAUSE_DECREMENT_SECONDS", "0.5") or "0.5"),
        adaptive_pause_max_seconds=float(get_env(env_values, "ERA5_ADAPTIVE_PAUSE_MAX_SECONDS", "30") or "30"),
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


def build_days(year: int, month: int) -> List[str]:
    last_day = calendar.monthrange(year, month)[1]
    return [f"{day:02d}" for day in range(1, last_day + 1)]


def build_target_path(output_dir: Path, variable: str, year: int, month: int, file_format: str) -> Path:
    ext = "grib" if file_format == "grib" else "nc"
    return output_dir / variable / f"era5_{variable}_{year}{month:02d}.{ext}"


def build_request(config: Config, variable: str, year: int, month: int) -> Dict[str, object]:
    return {
        "product_type": config.product_type,
        "variable": [variable],
        "year": f"{year:04d}",
        "month": f"{month:02d}",
        "day": build_days(year, month),
        "time": TIME_STEPS,
        "area": config.area,
        "data_format": config.file_format,
    }


def get_client(config: Config) -> cdsapi.Client:
    kwargs: Dict[str, object] = {"quiet": False, "progress": True}
    if config.cds_key:
        kwargs["url"] = config.cds_url or DEFAULT_CDS_API_URL
        kwargs["key"] = config.cds_key
    elif config.cds_url:
        kwargs["url"] = config.cds_url
    return cdsapi.Client(**kwargs)


def setup_logger(config: Config) -> logging.Logger:
    logger = logging.getLogger("era5_downloader")
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
    expected_suffix = ".grib" if file_format == "grib" else ".nc"
    if target.suffix.lower() != expected_suffix:
        return False, f"Estensione non coerente con file_format: {target.suffix} != {expected_suffix}"

    output_abs = output_dir.resolve()
    target_abs = target.resolve()
    try:
        target_abs.relative_to(output_abs)
    except ValueError:
        return False, f"Target fuori da ERA5_OUTPUT_DIR: {target_abs}"
    return True, "ok"


def validate_request_payload(
    config: Config,
    request: Dict[str, object],
    variable: str,
    year: int,
    month: int,
) -> Tuple[bool, str]:
    if request.get("product_type") != config.product_type:
        return False, "product_type non coerente"
    if request.get("data_format") != config.file_format:
        return False, "data_format non coerente"
    if request.get("year") != f"{year:04d}":
        return False, "year non coerente"
    if request.get("month") != f"{month:02d}":
        return False, "month non coerente"

    req_vars = request.get("variable")
    if not isinstance(req_vars, list) or req_vars != [variable]:
        return False, "variable non coerente"

    req_area = request.get("area")
    if not isinstance(req_area, list) or req_area != config.area:
        return False, "area non coerente"

    req_time = request.get("time")
    if not isinstance(req_time, list) or len(req_time) != 24:
        return False, "time non coerente (attese 24 ore)"

    req_days = request.get("day")
    if not isinstance(req_days, list) or len(req_days) != calendar.monthrange(year, month)[1]:
        return False, "day non coerente con il mese"

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
        if not (magic.startswith(b"CDF") or magic.startswith(b"\x89HDF")):
            return False, f"header NetCDF non valido: {magic!r}"

    return True, f"ok, size={size}"


def is_throttle_error(exc: Exception) -> bool:
    msg = str(exc).lower()
    throttle_patterns = [
        "429",
        "too many requests",
        "rate limit",
        "rate-limit",
        "throttl",
    ]
    return any(p in msg for p in throttle_patterns)


def download_single_chunk(
    client: cdsapi.Client,
    config: Config,
    logger: logging.Logger,
    variable: str,
    year: int,
    month: int,
) -> ChunkResult:
    ts = datetime.now().isoformat(timespec="seconds")
    target = build_target_path(config.output_dir, variable, year, month, config.file_format)
    target.parent.mkdir(parents=True, exist_ok=True)

    def _result(status: str, note: str, throttle_detected: bool = False) -> ChunkResult:
        size = target.stat().st_size if target.exists() else 0
        return ChunkResult(
            timestamp=ts,
            variable=variable,
            year=year,
            month=month,
            status=status,
            throttle_detected=throttle_detected,
            file_path=target,
            file_size_bytes=size,
            note=note,
        )

    target_ok, target_msg = validate_target_path(target, config.output_dir, config.file_format)
    if not target_ok:
        logger.error("[TARGET_INVALID] variable=%s year=%s month=%s target=%s msg=%s", variable, year, month, target, target_msg)
        return _result(STATUS_ERROR, f"TARGET_INVALID: {target_msg}")

    request = build_request(config, variable, year, month)
    req_ok, req_msg = validate_request_payload(config, request, variable, year, month)
    if not req_ok:
        logger.error("[REQUEST_INVALID] variable=%s year=%s month=%s msg=%s", variable, year, month, req_msg)
        return _result(STATUS_ERROR, f"REQUEST_INVALID: {req_msg}")

    logger.info(
        "[REQUEST_READY] dataset=%s variable=%s year=%s month=%02d format=%s area=%s target=%s",
        config.dataset,
        variable,
        year,
        month,
        config.file_format,
        config.area,
        target,
    )

    if config.skip_existing and target.exists() and target.stat().st_size > 0:
        ok, msg = verify_downloaded_file(target, config.file_format, config.verify_download_magic)
        if ok:
            logger.info("[SKIP_EXISTING_OK] variable=%s year=%s month=%02d target=%s %s", variable, year, month, target, msg)
            return _result(STATUS_SKIPPED, msg)
        logger.warning("[SKIP_EXISTING_BAD] variable=%s year=%s month=%02d target=%s msg=%s", variable, year, month, target, msg)
        target.unlink(missing_ok=True)

    if config.dry_run:
        logger.info("[DRY_RUN] variable=%s year=%s month=%02d target=%s", variable, year, month, target)
        return _result(STATUS_DRY_RUN, "dry_run")

    attempt = 0
    last_error = ""
    throttle_detected = False
    while attempt < config.max_retries:
        attempt += 1
        try:
            logger.info("[DOWNLOAD_START] variable=%s year=%s month=%02d attempt=%s/%s", variable, year, month, attempt, config.max_retries)
            client.retrieve(config.dataset, request, str(target))
            ok, msg = verify_downloaded_file(target, config.file_format, config.verify_download_magic)
            if not ok:
                logger.error("[VERIFY_FAILED] variable=%s year=%s month=%02d target=%s msg=%s", variable, year, month, target, msg)
                target.unlink(missing_ok=True)
                raise RuntimeError(msg)

            logger.info("[DOWNLOAD_DONE] variable=%s year=%s month=%02d target=%s %s", variable, year, month, target, msg)
            return _result(STATUS_DOWNLOADED, msg, throttle_detected=throttle_detected)
        except Exception as exc:  # cdsapi solleva eccezioni varie lato rete/server
            last_error = str(exc)
            current_is_throttle = is_throttle_error(exc)
            throttle_detected = throttle_detected or current_is_throttle

            logger.error("[DOWNLOAD_ERROR] variable=%s year=%s month=%02d attempt=%s/%s target=%s err=%s", variable, year, month, attempt, config.max_retries, target, exc)
            if current_is_throttle:
                logger.warning("[THROTTLE_DETECTED] variable=%s year=%s month=%02d attempt=%s/%s", variable, year, month, attempt, config.max_retries)

            if target.exists() and target.stat().st_size == 0:
                target.unlink(missing_ok=True)

            if attempt < config.max_retries:
                if current_is_throttle:
                    throttle_wait = min(
                        config.throttle_max_wait_seconds,
                        max(config.throttle_retry_wait_seconds, config.retry_wait_seconds) * (2 ** (attempt - 1)),
                    )
                    logger.info(
                        "[THROTTLE_WAIT] variable=%s year=%s month=%02d wait_seconds=%s",
                        variable,
                        year,
                        month,
                        throttle_wait,
                    )
                    time.sleep(throttle_wait)
                else:
                    logger.info("[RETRY_WAIT] variable=%s year=%s month=%02d wait_seconds=%s", variable, year, month, config.retry_wait_seconds)
                    time.sleep(config.retry_wait_seconds)

    return _result(STATUS_ERROR, f"max retries esauriti: {last_error}", throttle_detected=throttle_detected)


def _write_report(config: Config, results: List[ChunkResult], logger: logging.Logger) -> None:
    report_path = config.report_file
    report_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not report_path.exists()

    with report_path.open("a", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        if write_header:
            writer.writerow(["timestamp", "variable", "year", "month", "status", "throttle_detected", "file_path", "file_size_bytes", "note"])
        for r in results:
            writer.writerow([
                r.timestamp,
                r.variable,
                r.year,
                f"{r.month:02d}",
                r.status,
                str(r.throttle_detected).lower(),
                str(r.file_path),
                r.file_size_bytes,
                r.note,
            ])

    logger.info("[REPORT] Scritto: %s (%s righe)", report_path, len(results))


def run(config: Config) -> int:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(config)
    logger.info("=" * 80)
    logger.info("ERA5 hourly single levels downloader")
    logger.info("Dataset      : %s", config.dataset)
    logger.info("Periodo      : %s-%02d -> %s-%02d", config.start_year, config.start_month, config.end_year, config.end_month)
    logger.info("Area         : N/W/S/E = %s", config.area)
    logger.info("Output       : %s", config.output_dir)
    logger.info("Formato      : %s", config.file_format)
    logger.info("Variabili    : %s", len(config.variables))
    logger.info("Skip existing: %s", config.skip_existing)
    logger.info("Dry run      : %s", config.dry_run)
    logger.info("Log file     : %s", config.log_file)
    logger.info("Throttle wait: %ss (max %ss)", config.throttle_retry_wait_seconds, config.throttle_max_wait_seconds)
    logger.info("Adaptive pause: enabled=%s base=%ss max=%ss", config.adaptive_pause_enabled, config.request_pause_seconds, config.adaptive_pause_max_seconds)
    logger.info("=" * 80)

    client = get_client(config)

    results: List[ChunkResult] = []
    current_pause_seconds = config.request_pause_seconds
    for variable in config.variables:
        for year, month in month_iter(
            config.start_year,
            config.start_month,
            config.end_year,
            config.end_month,
        ):
            result = download_single_chunk(client, config, logger, variable, year, month)
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

            # Pausa solo se abbiamo effettivamente contattato il CDS (non per i chunk saltati).
            if result.status != STATUS_SKIPPED and current_pause_seconds > 0:
                logger.debug("[INTER_CHUNK_PAUSE] seconds=%s", round(current_pause_seconds, 3))
                time.sleep(current_pause_seconds)

    n_downloaded = sum(1 for r in results if r.status == STATUS_DOWNLOADED)
    n_skipped    = sum(1 for r in results if r.status == STATUS_SKIPPED)
    n_dry_run    = sum(1 for r in results if r.status == STATUS_DRY_RUN)
    n_error      = sum(1 for r in results if r.status == STATUS_ERROR)
    n_throttle   = sum(1 for r in results if r.throttle_detected)
    total        = len(results)

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
    parser = argparse.ArgumentParser(description="Download ERA5 hourly single levels per Europa")
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

    return run(config)


if __name__ == "__main__":
    raise SystemExit(main())
