#!/bin/bash -ue
if [ ! -f "/home/primeline/Data/IRnotator/test_out/manual_review/review.tsv" ]; then
    echo "ERROR: Manual review file not found:" >&2
    echo "  /home/primeline/Data/IRnotator/test_out/manual_review/review.tsv" >&2
    echo "" >&2
    echo "Export review.tsv from the topology reviewer (key e)," >&2
    echo "place it at that path, then re-run with -resume." >&2
    exit 1
fi
cp "/home/primeline/Data/IRnotator/test_out/manual_review/review.tsv" review.tsv
