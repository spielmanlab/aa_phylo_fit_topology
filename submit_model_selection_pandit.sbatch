#!/bin/bash
#SBATCH --job-name=modelselection
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-100


LINE=${SLURM_ARRAY_TASK_ID}


REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/dms_modelselection
ALNPATH=${TOPPATH}/pandit_aa_alignments
OUTPATH=${TOPPATH}/selected_models_pandit
mkdir -p $OUTPATH

FILE=${ALNPATH}/names.txt
INLINE=$(sed "${LINE}q;d" ${FILE})
IFS=',' read -r NAME <<< "$INLINE"

iqtree -nt 8 -m TESTONLY -s ${ALNPATH}/${NAME}.fasta -st AA
mv ${ALNPATH}/${NAME}.fasta.log $OUTPATH/${NAME}.model_selection_log
rm $ALNPATH/${NAME}.fasta.*


