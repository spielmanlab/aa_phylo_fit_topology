#!/bin/bash
#SBATCH --job-name=treeinf
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-10

REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/dms_modelselection/
ALNPATH=${TOPPATH}/simulations/alignments_empirical/
POGOPATH=${TOPPATH}/pogofit_models_empirical/
OUTPATH=${TOPPATH}/fitted_trees_empirical/
NAME=HA
cd $TOPPATH
for TREE in yeast greenplant greenalga opisthokonta prum ruhfel salichos dosreis anderson; do
    python3 infer_trees.py $NAME $TREE $REP $OUTPATH $ALNPATH $POGOPATH quantile_model_selection_empirical.csv 8
done


