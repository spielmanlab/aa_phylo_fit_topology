#!/bin/bash
#SBATCH --job-name=modelselection
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=11-20

NAME=HA
REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/dms_modelselection
ALNPATH=${TOPPATH}/simulations/alignments_empirical
OUTPATH=${TOPPATH}/selected_models_empirical
mkdir -p $OUTPATH
for TREE in yeast greenplant greenalga opisthokonta prum ruhfel salichos dosreis anderson; do
   X=${NAME}_${TREE}_rep${REP}_AA
   iqtree -nt 8 -m TESTONLY -s ${ALNPATH}/${X}.fasta -st AA
   mv ${ALNPATH}/${X}.fasta.log $OUTPATH/${X}.model_selection_log
   rm $ALNPATH/${X}.fasta.*
done

