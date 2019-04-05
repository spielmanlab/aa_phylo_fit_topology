#!/bin/bash
#SBATCH --job-name=pogomodel
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-10

REP=${SLURM_ARRAY_TASK_ID}
NAME=HA
source ~/.bash_profile

HYPATH=/csm_data/spielman_lab/bin/hyphy
LIBPATH=${HYPATH}/res
MPI=${HYPATH}/HYPHYMPI
MP=${HYPATH}/HYPHYMP
BF=${LIBPATH}/TemplateBatchFiles/pogofit_fixgamma.bf

OUTPATH=/csm_data/spielman_lab/dms_modelselection/pogofit_models_empirical
mkdir -p $OUTPATH

ALNPATH=/csm_data/spielman_lab/dms_modelselection/simulations/empirical_alignments
for TREE in yeast greenplant greenalga opisthokonta prum ruhfel salichos dosreis anderson; do
    DAT=${NAME}_${TREE}_rep${REP}_AA.dat
    echo $DAT
    if [ ! -f ${OUTPATH}/${DAT}.POGOFIT.json ]; then
         echo "INFERRING ${DAT}"
         $MP CPU=12 LIBPATH=$LIBPATH $BF Single ${ALNPATH}/$DAT JTT Yes PAML
         mv ${ALNPATH}/$DAT.* $OUTPATH
    fi
done
