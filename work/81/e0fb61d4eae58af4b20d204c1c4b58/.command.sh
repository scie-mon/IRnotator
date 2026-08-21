#!/bin/bash -ue
cd '/home/primeline/Data/IRnotator'
python '/home/primeline/Data/IRnotator/bin/validate_reference_config.py' \
    --project-dir '/home/primeline/Data/IRnotator' \
    --hmm-dir 'hmms' \
    --salvage-path '/home/primeline/Data/maker/metagenome/12_Dmagna_LRV'
