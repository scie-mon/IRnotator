#!/usr/bin/env python3
from pathlib import Path

source_path = Path("bin/irnotator.py")
target_path = Path("bin/irnotator.patched.py")
source = source_path.read_text(encoding="utf-8")

replacements = [
    (
        "from __future__ import annotations\n\nimport csv\n",
        "from __future__ import annotations\n\nimport ast\nimport csv\n",
    ),
    (
        "  irnotator -profile docker --proteins_faa candidates.faa\n",
        "  irnotator -profile docker --proteins_faa candidates.faa\n  irnotator --omit-warnings --outdir results\n",
    ),
    (
        "uses params.outdir resolved from `nextflow config -flat`.\n",
        "uses params.outdir resolved from `nextflow config -flat`.\n\n"
        "Use --omit-warnings to suppress reference-HMM and DeepTMHMM-cache\n"
        "configuration diagnostics.\n",
    ),
    (
        "command = [\"nextflow\", \"config\", \"-flat\", *nextflow_profile_options(arguments), str(PROJECT_DIR)]",
        "command = [\"nextflow\", \"config\", \"-flat\", *nextflow_config_options(arguments), str(PROJECT_DIR)]",
    ),
    (
        "    arguments = sys.argv[1:]\n    if arguments == [\"--help\"] or arguments == [\"-h\"]:\n",
        "    arguments = sys.argv[1:]\n"
        "    omit_warnings = \"--omit-warnings\" in arguments\n"
        "    arguments = [argument for argument in arguments if argument != \"--omit-warnings\"]\n"
        "    if arguments == [\"--help\"] or arguments == [\"-h\"]:\n",
    ),
    (
        "    launch_dir = Path.cwd().resolve()\n    outdir = parse_path_option(arguments, \"--outdir\", launch_dir) or resolve_config_outdir(arguments, launch_dir)\n",
        "    launch_dir = Path.cwd().resolve()\n"
        "    resolved_config = resolve_flat_config(arguments, launch_dir)\n"
        "    hmm_value = parse_value_option(arguments, \"--hmm_dir\") or unquote_config_value(resolved_config.get(\"params.hmm_dir\")) or \"hmms\"\n"
        "    salvage_value = parse_value_option(arguments, \"--deeptmhmm_salvage_paths\") or resolved_config.get(\"params.deeptmhmm_salvage_paths\") or \"[]\"\n"
        "    hmm_dir = resolve_runtime_path(hmm_value, launch_dir)\n"
        "    salvage_paths = parse_salvage_paths(salvage_value, launch_dir)\n"
        "    if not omit_warnings:\n"
        "        run_reference_audit(hmm_dir, salvage_paths, launch_dir)\n"
        "    outdir = parse_path_option(arguments, \"--outdir\", launch_dir) or resolve_config_outdir(arguments, launch_dir)\n",
    ),
    (
        "    command_base = [\"nextflow\", \"run\", str(PIPELINE), *arguments]\n",
        "    command_base = [\"nextflow\", \"run\", str(PIPELINE), *arguments, \"--omit_warnings\", \"true\"]\n",
    ),
]

helpers = r'''

def nextflow_config_options(arguments: list[str]) -> list[str]:
    options: list[str] = []
    index = 0
    options_with_values = {"-profile", "-c", "-C", "-params-file"}
    while index < len(arguments):
        argument = arguments[index]
        if argument in options_with_values:
            if index + 1 == len(arguments):
                fail(f"{argument} requires a value")
            options.extend([argument, arguments[index + 1]])
            index += 2
            continue
        if any(argument.startswith(f"{option}=") for option in options_with_values):
            option, value = argument.split("=", 1)
            if not value:
                fail(f"{option} requires a value")
            options.extend([option, value])
        index += 1
    return options


def resolve_flat_config(arguments: list[str], cwd: Path) -> dict[str, str]:
    command = ["nextflow", "config", "-flat", *nextflow_config_options(arguments), str(PROJECT_DIR)]
    result = subprocess.run(command, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip() or f"exit status {result.returncode}"
        fail(f"could not resolve pipeline configuration with {' '.join(command)}: {detail}")
    resolved: dict[str, str] = {}
    for line in result.stdout.splitlines():
        key, separator, value = line.partition("=")
        if separator:
            resolved[key.strip()] = value.strip()
    return resolved


def parse_value_option(arguments: list[str], option: str) -> str | None:
    value: str | None = None
    index = 0
    while index < len(arguments):
        argument = arguments[index]
        if argument == option:
            if index + 1 == len(arguments):
                fail(f"{option} requires a value")
            value = arguments[index + 1]
            index += 2
            continue
        prefix = f"{option}="
        if argument.startswith(prefix):
            value = argument[len(prefix):]
            if not value:
                fail(f"{option} requires a value")
        index += 1
    return value


def unquote_config_value(value: str | None) -> str | None:
    if value is None:
        return None
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        value = value[1:-1]
    return None if not value or value == "null" else value


def resolve_runtime_path(value: str, cwd: Path) -> Path:
    path = Path(unquote_config_value(value) or "").expanduser()
    return (cwd / path).resolve() if not path.is_absolute() else path.resolve()


def parse_salvage_paths(value: str, cwd: Path) -> list[Path]:
    value = (value or "").strip()
    if not value or value == "[]" or value == "null":
        return []
    try:
        parsed = ast.literal_eval(value)
    except (SyntaxError, ValueError):
        parsed = [item.strip().strip("'\"") for item in value.strip("[]").split(",") if item.strip()]
    if isinstance(parsed, str):
        parsed = [parsed]
    if not isinstance(parsed, (list, tuple)):
        fail(f"could not parse deeptmhmm_salvage_paths: {value}")
    return [resolve_runtime_path(str(path), cwd) for path in parsed]


def run_reference_audit(hmm_dir: Path, salvage_paths: list[Path], cwd: Path) -> None:
    command = [
        sys.executable,
        str(PROJECT_DIR / "bin" / "validate_reference_config.py"),
        "--project-dir", str(PROJECT_DIR),
        "--hmm-dir", str(hmm_dir),
    ]
    for path in salvage_paths:
        command.extend(["--salvage-path", str(path)])
    subprocess.run(command, cwd=cwd, check=False)
'''

anchor = "\ndef\ndef resolve_config_outdir(arguments: list[str], cwd: Path) -> Path:\n"
if anchor not in source:
    raise SystemExit("Could not locate the controller configuration helper anchor.")
source = source.replace(anchor, helpers + anchor, 1)
for old, new in replacements:
    if old not in source:
        raise SystemExit(f"Expected source fragment not found:\n{old}")
    source = source.replace(old, new, 1)

target_path.write_text(source, encoding="utf-8")
print(f"Wrote {target_path}")
