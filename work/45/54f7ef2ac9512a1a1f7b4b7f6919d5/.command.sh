#!/bin/bash -ue
awk '!/^#/ && NF > 0 {print $1}' PF00060.tbl | sort -u > HMM_ids.txt
