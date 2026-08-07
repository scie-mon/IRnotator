Manual review package
=====================

Identity rules
--------------
- seq_id is the IRnotator internal ID and is the only decision key.
- source_id/source_header are displayed only for human interpretation.
- image filenames use internal IDs.
- review.tsv must retain the exported seq_id values exactly.

Instructions
------------
1. From this directory, serve locally:
   python3 -m http.server 8000

2. Open:
   http://localhost:8000/topology-reviewer.html

3. Review every plot and export review.tsv with e.

4. Resume the pipeline after review.tsv has been written here.
