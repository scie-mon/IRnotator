Manual review package
=====================

The manifest is ordered by merged gene model, then isoform.

1. Serve this directory: python3 -m http.server 8000
2. Open topology-reviewer.html.
3. Export review.tsv into this directory.
4. Run ./rerun_after_review.sh /path/to/IRnotator/work
5. Re-run the same Nextflow command with -resume.

seq_id is the IRnotator internal ID and is the only decision key.
