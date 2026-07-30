process COLLECT_DEEPTMHMM_SUMMARY {
    tag "deeptmhmm_summary"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    val deeptmhmm_results

    output:
    path "deeptmhmm_summary", emit: summary_dir

    script:
    """
    mkdir -p deeptmhmm_summary

    python - <<'PY'
from pathlib import Path

raw = ${deeptmhmm_results.collect { "'${it}'" }}

# Collected list of (seq_id, path) tuples flattens to [id, path, id, path, ...]
paths = []
i = 0
while i < len(raw):
    p = Path(str(raw[i]))
    if p.exists() and p.is_dir():
        paths.append(p)
        i += 1
        continue
    if i + 1 < len(raw):
        p2 = Path(str(raw[i + 1]))
        if not p2.exists():
            raise FileNotFoundError(f"DeepTMHMM result path does not exist: {p2}")
        paths.append(p2)
        i += 2
        continue
    raise ValueError(f"Unpaired item in deeptmhmm_results collect: {raw[i]!r}")

if not paths:
    raise RuntimeError("No DeepTMHMM result directories to merge")

out = Path("deeptmhmm_summary/TMRs.gff3")
wrote_gff_version = False
with out.open("w") as w:
    for idx, d in enumerate(paths):
        gff = d / "TMRs.gff3"
        if not gff.is_file():
            raise FileNotFoundError(f"missing TMRs.gff3 in {d}")

        # Batch-mode DeepTMHMM separates sequence blocks with //
        if idx > 0:
            w.write("//")
            w.write(chr(10))

        with gff.open() as fh:
            for line in fh:
                if line.startswith("##gff-version"):
                    if not wrote_gff_version:
                        w.write(line)
                        wrote_gff_version = True
                    continue
                # Keep per-sequence # headers, feature lines, and any // already present
                w.write(line)
PY
    """
}
