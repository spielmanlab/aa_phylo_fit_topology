#!/bin/bash
#SBATCH --job-name=sh_sim
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-20

REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/dms_modelselection
SCRIPTPATH=$TOPPATH/run_analysis_files

cd $SCRIPTPATH
python3 run_sh_tests_poisson.py simulation $REP

