#!/bin/bash
#SBATCH --job-name=pogomodel
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-100

LINE=${SLURM_ARRAY_TASK_ID}

REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/dms_modelselection
ALNPATH=${TOPPATH}/pandit_aa_alignments
OUTPATH=${TOPPATH}/pogofit_models_pandit
mkdir -p $OUTPATH

FILE=${ALNPATH}/names.txt
INLINE=$(sed "${LINE}q;d" ${FILE})
IFS=',' read -r NAME <<< "$INLINE"


HYPATH=/csm_data/spielman_lab/bin/hyphy
LIBPATH=${HYPATH}/res
MPI=${HYPATH}/HYPHYMPI
MP=${HYPATH}/HYPHYMP
BF=${LIBPATH}/TemplateBatchFiles/pogofit_fixgamma.bf


DAT=${NAME}.dat
echo $DAT
if [ ! -f ${OUTPATH}/${DAT}.POGOFIT.json ]; then
    echo "INFERRING ${DAT}"
    $MP CPU=12 LIBPATH=$LIBPATH $BF Single ${ALNPATH}/$DAT JTT Yes PAML
    mv ${ALNPATH}/$DAT.* $OUTPATH
fi
