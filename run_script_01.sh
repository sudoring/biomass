#!/bin/bash

module load python3/3.10.12

echo Submitting $1
python -u  /usr3/graduate/mkmoon/GitHub/biomass/download_script_multi.py $1
