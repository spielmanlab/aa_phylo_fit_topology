#!/bin/bash
#SBATCH --job-name=treeinf
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-20

REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/dms_modelselection/
ALNPATH=${TOPPATH}/simulations/alignments_empirical/
POGOPATH=${TOPPATH}/pogofit_models_empirical/
OUTPATH=${TOPPATH}/fitted_trees_empirical/
NAME=$1
cd $TOPPATH
for TREE in greenalga opisthokonta prum ruhfel salichos dosreis andersen rayfinned spiralia; do
   python3 infer_trees.py $NAME $TREE $REP $OUTPATH $ALNPATH $POGOPATH processed_model_selection/quantile_model_selection_empirical.csv 8
done


