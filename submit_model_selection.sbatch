#!/bin/bash
#SBATCH --job-name=modelselection
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-20

REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/dms_modelselection
ALNPATH=${TOPPATH}/simulations/alignments
OUTPATH=${TOPPATH}/selected_models
mkdir -p $OUTPATH
for TREE in rtree100_bl0.03_rep1; do ####rtree100_bl0.3_rep1 rtree100_bl1.5_rep1 rtree100_bl3_rep1; do
  for NAME in HIV LAC Gal4 NP HA; do
   X=${NAME}_${TREE}_rep${REP}_AA
   iqtree -nt 8 -m TESTONLY -s ${ALNPATH}/${X}.fasta -st AA
   mv ${ALNPATH}/${X}.fasta.log $OUTPATH/${X}.model_selection_log
   rm $ALNPATH/${X}.fasta.*
  done
done

