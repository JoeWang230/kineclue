#!/bin/bash
# Full run: analytical stage -> fill five-frequency values -> numerical stage
set -e
cd "$(dirname "$0")"
K=../..
python3 $K/kineclue_main.py main_AlSn_input.txt
python3 make_5fm.py CALC/
python3 $K/kineclue_num.py num_AlSn_input.txt
python3 $K/kineclue_num.py num_AlSn_lowT_input.txt  # 300-800 K, for D_v in Al-0.04 at.% Sn
