#!/bin/bash
#SBATCH --job-name=sh_pandit
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-200

LINE=${SLURM_ARRAY_TASK_ID}

REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/dms_modelselection
SCRIPTPATH=$TOPPATH/run_analysis_files
ALNPATH=$TOPPATH/pandit_aa_alignments/
INFPATH=$TOPPATH/fitted_trees_pandit/
QFILE=$TOPPATH/processed_model_selection/quantile_model_selection_pandit.csv


FILE=${ALNPATH}/names.txt
INLINE=$(sed "${LINE}q;d" ${FILE})
IFS=',' read -r NAME <<< "$INLINE"



cd $SCRIPTPATH
python3 run_sh_tests.py pandit $NAME $ALNPATH $INFPATH $QFILE
