#!/bin/bash
#SBATCH --job-name=modelselection
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-5

REP=${SLURM_ARRAY_TASK_ID}
NAME=$1
TOPPATH=/csm_data/spielman_lab/dms_modelselection/
ALNPATH=${TOPPATH}/simulations/alignments
OUTPATH=${TOPPATH}/selected_models
mkdir -p $OUTPATH
for TREE in rtree64_bl0.3.tree rtree64_bl3.0.tree btree64_bl0.3.tree btree64_bl3.0.tree; do
    iqtree -nt 8 -m TESTONLY -s $ALNPATH/${NAME}_${TREE}_rep${REP}_AA.fasta -st AA
    mv $ALNPATH/${NAME}_${TREE}_rep${REP}_AA.fasta.log $OUTPATH/${NAME}_${TREE}_rep${REP}_AA.model_selection_log
    rm $ALNPATH/${NAME}_${TREE}_rep${REP}_AA.fasta.*
done
#HA_rtree64_bl0.3.tree_rep4_AA.fasta   


