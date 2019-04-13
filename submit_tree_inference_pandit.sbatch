#!/bin/bash
#SBATCH --job-name=treeinf_pandit
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-100

LINE=${SLURM_ARRAY_TASK_ID}

REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/dms_modelselection
ALNPATH=${TOPPATH}/pandit_aa_alignments/
POGOPATH=${TOPPATH}/pogofit_models_pandit/
OUTPATH=${TOPPATH}/fitted_trees_pandit/
mkdir -p $OUTPATH

FILE=${ALNPATH}/names.txt
INLINE=$(sed "${LINE}q;d" ${FILE})
IFS=',' read -r NAME <<< "$INLINE"

cd $TOPPATH
python3 infer_trees.py $NAME NA NA $OUTPATH $ALNPATH $POGOPATH processed_model_selection/quantile_model_selection_pandit.csv 8
