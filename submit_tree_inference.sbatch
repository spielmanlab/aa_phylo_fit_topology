#!/bin/bash
#SBATCH --job-name=modelselection
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-20

REP=${SLURM_ARRAY_TASK_ID}
NAME=$1
TOPPATH=/csm_data/spielman_lab/dms_modelselection/
ALNPATH=${TOPPATH}/simulations/alignments/
POGOPATH=${TOPPATH}/pogofit_models/
OUTPATH=${TOPPATH}/fitted_trees/
mkdir -p $OUTPATH

cd $TOPPATH
for TREE in rtree100_bl0.3_rep1 rtree100_bl0.75_rep1 rtree100_bl1.5_rep1 rtree100_bl3_rep1; do
    python3 infer_trees.py $NAME $TREE $REP $OUTPATH $ALNPATH $POGOPATH quantile_model_selection.csv 8
done
