#!/bin/bash -ue
if grep -q '^>' unresolved_deeptmhmm.faa; then
    echo 'DeepTMHMM results are missing for one or more sequences:' >&2
    cat unresolved_deeptmhmm.faa >&2
    exit 1
fi
