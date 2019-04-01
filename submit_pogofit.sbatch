#!/bin/bash
#SBATCH --job-name=pogomodel
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-5

REP=${SLURM_ARRAY_TASK_ID}

source ~/.bash_profile

HYPATH=/csm_data/spielman_lab/bin/hyphy
LIBPATH=${HYPATH}/res
MPI=${HYPATH}/HYPHYMPI
MP=${HYPATH}/HYPHYMP
BF=${LIBPATH}/TemplateBatchFiles/pogofit_fixgamma.bf

OUTPATH=/csm_data/spielman_lab/dms_modelselection/pogofit_models
mkdir -p $OUTPATH

ALNPATH=/csm_data/spielman_lab/dms_modelselection/simulations/alignments
NAME=$1
#PF00135 PF00171 PF07690 PF00496 PF01094
#for NAME in NP LAC HA HIV Gal4; do
  for TREE in rtree64_bl0.3.tree rtree64_bl3.0.tree btree64_bl0.3.tree btree64_bl3.0.tree; do
    DAT=${NAME}_${TREE}_rep${REP}_AA.dat
    echo $DAT
    if [ ! -f ${OUTPATH}/${DAT}.POGOFIT.json ]; then
         echo "INFERRING ${DAT}"
         $MP CPU=12 LIBPATH=$LIBPATH $BF Single ${ALNPATH}/$DAT JTT Yes PAML
         mv ${ALNPATH}/$DAT.* $OUTPATH
    fi
  done
#done
