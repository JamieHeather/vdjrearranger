from datetime import datetime
import time
from pathlib import Path
import vdjrearranger


def write_run_log(outdir: Path, parameters: dict, start_time: float):
    """
    Writes execution metadata, timestamp, runtime, and CLI parameters to a log file.

    :param outdir: str, directory path where the run.log will be saved.
    :param parameters: dict containing all execution flags and their states.
    :param start_time: float timestamp generated when execution began.
    """
    log_path = outdir / "run.log"
    now = datetime.now()
    duration = time.time() - start_time

    with open(log_path, "w") as handle:
        handle.write(f"vdjrearranger_version: {getattr(vdjrearranger, '__version__', 'unknown')}\n")
        handle.write(f"run_date: {now.strftime('%Y-%m-%d')}\n")
        handle.write(f"run_time: {now.strftime('%H:%M:%S')}\n")
        handle.write(f"run_duration_seconds: {duration:.2f}\n\n")
        handle.write("=== Parameters ===\n")
        for key, value in parameters.items():
            handle.write(f"{key}: {value}\n")
