#!/bin/bash
#SBATCH --job-name=treeinf
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-20

REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/dms_modelselection/results_analysis

cd $TOPPATH
python3 run_sh_tests.py empirical $REP
