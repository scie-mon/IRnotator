Manual review package
=====================

1. From this directory, serve locally:
   python3 -m http.server 8000

2. Open:
   http://localhost:8000/topology-reviewer.html

3. Review every plot. Export review.tsv with key e.

4. Resume later with:
   --review_tsv /path/to/review.tsv -resume
