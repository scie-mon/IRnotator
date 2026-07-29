#!/bin/bash -ue
hmmsearch         -E 1e-5         --tblout PF10613.tbl         PF10613.hmm         test.faa         > PF10613.txt
