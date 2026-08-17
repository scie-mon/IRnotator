#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import signal
import socket
import subprocess
import sys
import threading
import time
import webbrowser
from datetime import datetime, timezone
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import quote

SCRIPT_PATH = Path(__file__).resolve()
PROJECT_DIR = SCRIPT_PATH.parent.parent
PIPELINE = PROJECT_DIR / "main.nf"
SESSION_NAME = ".irnotator-review-session.json"


def fail(message: str) -> None:
    raise SystemExit(f"irnotator: ERROR: {message}")


def usage() -> str:
    return """Usage:
  irnotator [NEXTFLOW_OPTIONS_AND_PARAMETERS]

Runs main.nf once, serves the generated manual-review package locally, then
runs successive review rounds after you export and copy review.tsv.

Examples:
  irnotator -profile docker --proteins_faa candidates.faa --outdir results
  irnotator -profile docker --genome_fasta genome.fa --annot_gff 'input/*.gff3' --outdir results

All arguments are passed unchanged to `nextflow run main.nf`, except -resume,
which is managed by this controller.
"""


def parse_path_option(arguments: list[str], option: str, default: Path, cwd: Path) -> Path:
    for index, argument in enumerate(arguments):
        if argument == option:
            if index + 1 == len(arguments):
                fail(f"{option} requires a value")
            value = arguments[index + 1]
            path = Path(value).expanduser()
            return (cwd / path).resolve() if not path.is_absolute() else path.resolve()
        prefix = f"{option}="
        if argument.startswith(prefix):
            value = argument[len(prefix):]
            if not value:
                fail(f"{option} requires a value")
            path = Path(value).expanduser()
            return (cwd / path).resolve() if not path.is_absolute() else path.resolve()
    return default


def find_free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("127.0.0.1", 0))
        return int(probe.getsockname()[1])


def write_session(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n")


def read_manifest(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"cannot read reviewer manifest {path}: {exc}")


def validate_review(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"review TSV is missing or empty: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            fail(f"review TSV has no header: {path}")
        required = {"decision"}
        if not required.issubset(reader.fieldnames) or not ({"seq_id", "protein_id"} & set(reader.fieldnames)):
            fail(f"review TSV requires decision and seq_id (or protein_id) columns: {path}")
        for line_number, row in enumerate(reader, 2):
            sequence_id = (row.get("seq_id") or row.get("protein_id") or "").strip()
            decision = (row.get("decision") or "").strip()
            if not sequence_id:
                fail(f"review TSV line {line_number} has no sequence ID")
            if decision not in {"accept", "reject"}:
                fail(f"review TSV line {line_number} has invalid decision {decision!r}")


class NoCacheHandler(SimpleHTTPRequestHandler):
    def end_headers(self) -> None:
        self.send_header("Cache-Control", "no-store, max-age=0")
        super().end_headers()

    def log_message(self, format: str, *args: object) -> None:
        return


def start_server(directory: Path, port: int) -> tuple[ThreadingHTTPServer, threading.Thread]:
    handler = lambda *args, **kwargs: NoCacheHandler(*args, directory=str(directory), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", port), handler)
    thread = threading.Thread(target=server.serve_forever, name="irnotator-review-server", daemon=True)
    thread.start()
    return server, thread


def open_reviewer(port: int, revision: object) -> str:
    url = (
        f"http://127.0.0.1:{port}/topology-reviewer.html"
        f"?revision={quote(str(revision))}&t={time.time_ns()}"
    )

    print(f"Reviewer opened for revision {revision}.")
    print(f"  {url}")

    try:
        subprocess.Popen(
            ["xdg-open", url],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
        )
    except FileNotFoundError:
        if not webbrowser.open(url):
            print("Could not open a browser automatically; open the URL above manually.", file=sys.stderr)

    return url


def run_nextflow(command: list[str], cwd: Path, log_path: Path) -> None:
    print("\nRunning:")
    print("  " + " ".join(command))
    print()

    with log_path.open("a") as log:
        log.write(
            f"\n[{datetime.now(timezone.utc).isoformat()}] "
            f"{' '.join(command)}\n"
        )

    result = subprocess.run(command, cwd=cwd)

    with log_path.open("a") as log:
        log.write(
            f"[{datetime.now(timezone.utc).isoformat()}] "
            f"exit status: {result.returncode}\n"
        )

    if result.returncode != 0:
        fail(f"Nextflow exited with status {result.returncode}; see {log_path}")


def prompt_for_review(review_path: Path, manifest_path: Path, expected_revision: object) -> None:
    manifest_mtime = manifest_path.stat().st_mtime_ns
    while True:
        print("\nExport review.tsv from the browser and copy it to:")
        print(f"  {review_path}")
        command = input("Type 'continue' when ready, 'open' to reopen the reviewer, or 'quit': ").strip().lower()
        if command == "quit":
            raise KeyboardInterrupt
        if command == "open":
            continue
        if command != "continue":
            print("Expected: continue, open, or quit.")
            continue
        if not review_path.exists():
            print("review.tsv is not present yet.")
            continue
        if review_path.stat().st_mtime_ns < manifest_mtime:
            print("review.tsv predates the current review package; export a fresh file.")
            continue
        try:
            validate_review(review_path)
        except SystemExit as exc:
            print(exc, file=sys.stderr)
            continue
        return


def main() -> None:
    arguments = sys.argv[1:]
    if not arguments or arguments == ["--help"] or arguments == ["-h"]:
        print(usage())
        return
    if any(argument == "-resume" or argument.startswith("-resume=") for argument in arguments):
        fail("do not supply -resume; irnotator controls review-round resumption")
    if not PIPELINE.is_file():
        fail(f"pipeline script not found: {PIPELINE}")

    launch_dir = Path.cwd().resolve()
    outdir = parse_path_option(arguments, "--outdir", launch_dir / "results", launch_dir)
    workdir = parse_path_option(arguments, "-work-dir", launch_dir / "work", launch_dir)
    review_dir = outdir / "manual_review"
    review_path = review_dir / "review.tsv"
    manifest_path = review_dir / "review-manifest.json"
    helper_path = review_dir / "rerun_after_review.sh"
    session_path = outdir / SESSION_NAME
    log_path = outdir / "irnotator-review.log"

    if review_path.exists():
        fail(f"existing review TSV found: {review_path}. Remove it or apply it manually before starting a new controller session.")

    outdir.mkdir(parents=True, exist_ok=True)
    command_base = ["nextflow", "run", str(PIPELINE), *arguments]
    session = {
        "project_dir": str(PROJECT_DIR),
        "launch_dir": str(launch_dir),
        "outdir": str(outdir),
        "workdir": str(workdir),
        "nextflow_command": command_base,
        "started_at": datetime.now(timezone.utc).isoformat(),
        "rounds_completed": 0,
    }
    write_session(session_path, session)

    run_nextflow(command_base, launch_dir, log_path)
    if not manifest_path.is_file() or not helper_path.is_file():
        fail(f"pipeline completed but no usable manual-review package was published in {review_dir}")

    port = find_free_port()
    server, _ = start_server(review_dir, port)
    print(f"Serving {review_dir} at http://127.0.0.1:{port}/")

    try:
        while True:
            manifest = read_manifest(manifest_path)
            revision = manifest.get("reviewer_revision", "unknown")
            session["current_reviewer_revision"] = revision
            session["current_manifest"] = str(manifest_path)
            write_session(session_path, session)
            open_reviewer(port, revision)
            prompt_for_review(review_path, manifest_path, revision)

            print(f"Applying {review_path}.")
            subprocess.run([str(helper_path), str(workdir)], cwd=review_dir, check=True)
            run_nextflow([*command_base, "-resume"], launch_dir, log_path)

            session["rounds_completed"] += 1
            session["last_round_completed_at"] = datetime.now(timezone.utc).isoformat()
            write_session(session_path, session)

            if review_path.exists():
                fail(f"the resumed pipeline left a root review TSV in place: {review_path}")
            if not manifest_path.is_file():
                fail(f"resumed pipeline did not publish {manifest_path}")
            print("\nNext review package is ready.")
    except KeyboardInterrupt:
        print("\nController stopped. The Nextflow outputs and review package were preserved.")
    finally:
        server.shutdown()
        server.server_close()


if __name__ == "__main__":
    main()
