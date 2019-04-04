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
ALNPATH=${TOPPATH}/simulations/alignments
OUTPATH=${TOPPATH}/selected_models
mkdir -p $OUTPATH

cd $TOPPATH
for TREE in rtree100_bl0.3_rep1 rtree100_bl1.5_rep1 rtree100_bl3_rep1; do
    python3 infer_trees.py $NAME $TREE $REP 8
done
#HA_rtree64_bl0.3.tree_rep4_AA.fasta   


