#!/bin/bash
#SBATCH --job-name=treeinf
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-20

REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/dms_modelselection
SCRIPTPATH=$TOPPATH/run_analysis_files
ALNPATH=$TOPPATH/simulations/alignments/
INFPATH=$TOPPATH/fitted_trees_simulation/
QFILE=$TOPPATH/processed_model_selection/quantile_model_selection_simulation.csv

cd $SCRIPTPATH
python3 run_sh_tests.py simulation $REP $ALNPATH $INFPATH $QFILE

