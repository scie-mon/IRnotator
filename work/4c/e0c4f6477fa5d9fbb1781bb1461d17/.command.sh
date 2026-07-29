#!/bin/bash -ue
export BIOLIB_DOCKER_RUNTIME=nvidia && biolib run --local 'DTU/DeepTMHMM:1.0.24' --fasta HMM_hits.aa
