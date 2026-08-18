#!/usr/bin/env python3
from __future__ import annotations

import csv
import io
import json
import os
import re
import secrets
import select
import socket
import subprocess
import sys
import tempfile
import threading
import time
import webbrowser
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import quote, urlparse

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
runs successive review rounds after you export, copy, or submit review.tsv.

Docker is the default profile. Supply -profile explicitly to override it.

Examples:
  irnotator
  irnotator --outdir results
  irnotator -profile apptainer
  irnotator -profile docker --proteins_faa candidates.faa

All arguments are passed unchanged to `nextflow run main.nf`, except -resume,
which is managed by this controller. When --outdir is omitted, the controller
uses params.outdir resolved from `nextflow config -flat`.
"""


def parse_path_option(arguments: list[str], option: str, cwd: Path) -> Path | None:
    for index, argument in enumerate(arguments):
        if argument == option:
            if index + 1 == len(arguments):
                fail(f"{option} requires a value")
            value = arguments[index + 1]
        else:
            prefix = f"{option}="
            if not argument.startswith(prefix):
                continue
            value = argument[len(prefix):]
            if not value:
                fail(f"{option} requires a value")
        path = Path(value).expanduser()
        return (cwd / path).resolve() if not path.is_absolute() else path.resolve()
    return None


def nextflow_profile_options(arguments: list[str]) -> list[str]:
    options: list[str] = []
    index = 0
    while index < len(arguments):
        argument = arguments[index]
        if argument == "-profile":
            if index + 1 == len(arguments):
                fail("-profile requires a value")
            options.extend(["-profile", arguments[index + 1]])
            index += 2
            continue
        if argument.startswith("-profile="):
            profile = argument.removeprefix("-profile=")
            if not profile:
                fail("-profile requires a value")
            options.extend(["-profile", profile])
        index += 1
    return options


def resolve_config_outdir(arguments: list[str], cwd: Path) -> Path:
    command = ["nextflow", "config", "-flat", *nextflow_profile_options(arguments), str(PROJECT_DIR)]
    result = subprocess.run(command, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip() or f"exit status {result.returncode}"
        fail(f"could not resolve params.outdir with {' '.join(command)}: {detail}")
    match = re.search(r"^params\.outdir\s*=\s*(.+?)\s*$", result.stdout, flags=re.MULTILINE)
    if not match:
        fail("--outdir was not supplied and `nextflow config -flat` did not resolve params.outdir; pass --outdir explicitly")
    value = match.group(1).strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        value = value[1:-1]
    if not value or value == "null":
        fail("--outdir was not supplied and resolved params.outdir is empty or null; pass --outdir explicitly")
    path = Path(value).expanduser()
    return (cwd / path).resolve() if not path.is_absolute() else path.resolve()


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


def validate_review_text(text: str, label: str) -> None:
    if not text.strip():
        fail(f"review TSV is missing or empty: {label}")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if not reader.fieldnames:
        fail(f"review TSV has no header: {label}")
    if "decision" not in reader.fieldnames or not ({"seq_id", "protein_id"} & set(reader.fieldnames)):
        fail(f"review TSV requires decision and seq_id (or protein_id) columns: {label}")
    for line_number, row in enumerate(reader, 2):
        sequence_id = (row.get("seq_id") or row.get("protein_id") or "").strip()
        decision = (row.get("decision") or "").strip()
        if not sequence_id:
            fail(f"review TSV line {line_number} has no sequence ID")
        if decision not in {"accept", "reject"}:
            fail(f"review TSV line {line_number} has invalid decision {decision!r}")


def validate_review(path: Path) -> None:
    if not path.is_file():
        fail(f"review TSV is missing: {path}")
    validate_review_text(path.read_text(), str(path))


class ReviewServer(ThreadingHTTPServer):
    def __init__(self, address: tuple[str, int], handler, directory: Path, review_path: Path):
        super().__init__(address, handler)
        self.directory = directory
        self.review_path = review_path
        self.reviewer_revision: object | None = None
        self.token = ""
        self.submitted = threading.Event()
        self.lock = threading.Lock()

    def begin_round(self, revision: object) -> None:
        with self.lock:
            self.reviewer_revision = revision
            self.token = secrets.token_urlsafe(32)
            self.submitted.clear()


class ReviewHandler(SimpleHTTPRequestHandler):
    server: ReviewServer

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(kwargs.pop("directory")), **kwargs)

    def end_headers(self) -> None:
        self.send_header("Cache-Control", "no-store, max-age=0")
        super().end_headers()

    def log_message(self, format: str, *args: object) -> None:
        return

    def send_json(self, status: HTTPStatus, payload: dict) -> None:
        body = json.dumps(payload).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        if urlparse(self.path).path == "/reviewer-session.json":
            with self.server.lock:
                payload = {
                    "token": self.server.token,
                    "reviewer_revision": self.server.reviewer_revision,
                }
            self.send_json(HTTPStatus.OK, payload)
            return
        super().do_GET()

    def do_POST(self) -> None:
        if urlparse(self.path).path != "/submit-review":
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        try:
            content_length = int(self.headers.get("Content-Length", "0"))
            if content_length <= 0 or content_length > 20_000_000:
                raise ValueError("invalid request size")
            payload = json.loads(self.rfile.read(content_length))
            if not isinstance(payload, dict):
                raise ValueError("request must be a JSON object")
            token = payload.get("token")
            revision = payload.get("reviewer_revision")
            tsv = payload.get("tsv")
            if not isinstance(token, str) or not isinstance(tsv, str):
                raise ValueError("request requires token and tsv strings")
            with self.server.lock:
                if token != self.server.token:
                    raise ValueError("invalid reviewer session token")
                if revision != self.server.reviewer_revision:
                    raise ValueError("reviewer revision does not match the active package")
                if self.server.submitted.is_set() or self.server.review_path.exists():
                    raise ValueError("a review TSV has already been submitted for this round")
                validate_review_text(tsv, "submitted review TSV")
                self.server.review_path.parent.mkdir(parents=True, exist_ok=True)
                with tempfile.NamedTemporaryFile(
                    mode="w",
                    encoding="utf-8",
                    newline="",
                    dir=self.server.review_path.parent,
                    prefix=".review.submit.",
                    suffix=".tsv",
                    delete=False,
                ) as handle:
                    handle.write(tsv)
                    temp_path = Path(handle.name)
                os.replace(temp_path, self.server.review_path)
                self.server.submitted.set()
        except (ValueError, json.JSONDecodeError, OSError, SystemExit) as exc:
            self.send_json(HTTPStatus.BAD_REQUEST, {"ok": False, "error": str(exc)})
            return
        self.send_json(HTTPStatus.OK, {"ok": True, "message": "Review submitted. Starting the next pipeline round."})


def start_server(directory: Path, review_path: Path, port: int) -> ReviewServer:
    handler = lambda *args, **kwargs: ReviewHandler(*args, directory=directory, **kwargs)
    server = ReviewServer(("127.0.0.1", port), handler, directory, review_path)
    threading.Thread(target=server.serve_forever, name="irnotator-review-server", daemon=True).start()
    return server


def open_reviewer(port: int, revision: object) -> str:
    url = f"http://127.0.0.1:{port}/topology-reviewer.html?revision={quote(str(revision))}&t={time.time_ns()}"
    print(f"Reviewer opened for revision {revision}.")
    print(f"  {url}")
    try:
        subprocess.Popen(["xdg-open", url], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
    except FileNotFoundError:
        if not webbrowser.open(url):
            print("Could not open a browser automatically; open the URL above manually.", file=sys.stderr)
    return url


def run_nextflow(command: list[str], cwd: Path, log_path: Path) -> None:
    print("\nRunning:")
    print("  " + " ".join(command))
    print()
    with log_path.open("a") as log:
        log.write(f"\n[{datetime.now(timezone.utc).isoformat()}] {' '.join(command)}\n")
    result = subprocess.run(command, cwd=cwd)
    with log_path.open("a") as log:
        log.write(f"[{datetime.now(timezone.utc).isoformat()}] exit status: {result.returncode}\n")
    if result.returncode != 0:
        fail(f"Nextflow exited with status {result.returncode}; see {log_path}")


def prompt_for_review(review_path: Path, server: ReviewServer, port: int, revision: object) -> None:
    print("\nSubmit from the browser, or export/copy review.tsv to:")
    print(f"  {review_path}")
    print("Type 'continue' to apply a copied TSV, 'open' to reopen the reviewer, or 'quit' to stop.")
    while True:
        if server.submitted.wait(timeout=0.25):
            print("Browser review submitted.")
            return
        readable, _, _ = select.select([sys.stdin], [], [], 0)
        if not readable:
            continue
        command = sys.stdin.readline()
        if command == "":
            raise KeyboardInterrupt
        command = command.strip().lower()
        if command == "quit":
            raise KeyboardInterrupt
        if command == "open":
            open_reviewer(port, revision)
            continue
        if command != "continue":
            print("Expected: continue, open, or quit.")
            continue
        if not review_path.exists():
            print("review.tsv is not present yet.")
            continue
        try:
            validate_review(review_path)
        except SystemExit as exc:
            print(exc, file=sys.stderr)
            continue
        return


def main() -> None:
    arguments = sys.argv[1:]
    if arguments == ["--help"] or arguments == ["-h"]:
        print(usage())
        return
    if not any(argument == "-profile" or argument.startswith("-profile=") for argument in arguments):
        arguments = ["-profile", "docker", *arguments]
    if any(argument == "-resume" or argument.startswith("-resume=") for argument in arguments):
        fail("do not supply -resume; irnotator controls review-round resumption")
    if not PIPELINE.is_file():
        fail(f"pipeline script not found: {PIPELINE}")

    launch_dir = Path.cwd().resolve()
    outdir = parse_path_option(arguments, "--outdir", launch_dir) or resolve_config_outdir(arguments, launch_dir)
    workdir = parse_path_option(arguments, "-work-dir", launch_dir) or (launch_dir / "work")
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
        "project_dir": str(PROJECT_DIR), "launch_dir": str(launch_dir), "outdir": str(outdir),
        "workdir": str(workdir), "nextflow_command": command_base,
        "started_at": datetime.now(timezone.utc).isoformat(), "rounds_completed": 0,
    }
    write_session(session_path, session)

    run_nextflow(command_base, launch_dir, log_path)
    if not manifest_path.is_file() or not helper_path.is_file():
        fail(f"pipeline completed but no usable manual-review package was published in {review_dir}")

    port = find_free_port()
    server = start_server(review_dir, review_path, port)
    print(f"Serving {review_dir} at http://127.0.0.1:{port}/")
    try:
        while True:
            manifest = read_manifest(manifest_path)
            revision = manifest.get("reviewer_revision", "unknown")
            server.begin_round(revision)
            session["current_reviewer_revision"] = revision
            session["current_manifest"] = str(manifest_path)
            write_session(session_path, session)
            open_reviewer(port, revision)
            prompt_for_review(review_path, server, port, revision)

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
