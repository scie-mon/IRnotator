process SPLIT_FASTA_BY_SEQID {
    tag "split_fasta_by_seqid"

    input:
    path fasta

    output:
    path "*.faa", emit: fasta_files

    script:
    """
    python - <<'PY'
    from pathlib import Path
    from Bio import SeqIO
    import re

    inp = Path("${fasta}")
    for rec in SeqIO.parse(str(inp), "fasta"):
        seq_id = re.sub(r'[^A-Za-z0-9_.-]+', '_', rec.id)
        out = Path(f"{seq_id}.faa")
        SeqIO.write(rec, str(out), "fasta")
    PY
    """
}
