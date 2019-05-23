#!/bin/bash
#SBATCH --job-name=sh_pandit
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-200

LINE=${SLURM_ARRAY_TASK_ID}

REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/dms_modelselection/results_analysis
mkdir -p $OUTPATH

FILE=${ALNPATH}/names.txt
INLINE=$(sed "${LINE}q;d" ${FILE})
IFS=',' read -r NAME <<< "$INLINE"

cd $TOPPATH
python3 run_sh_tests.py pandit $NAME