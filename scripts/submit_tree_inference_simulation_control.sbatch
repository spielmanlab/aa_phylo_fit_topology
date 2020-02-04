#!/bin/bash
#SBATCH --job-name=treeinf
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-20

REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/aa_phylo_fit_topology
ALNPATH=${TOPPATH}/simulations/alignments/wag_control/
OUTPATH=${TOPPATH}/fitted_trees_ufb_simulation_control/
mkdir -p ${OUTPATH}
NAME=$1
cd $TOPPATH
for TREE in opisthokonta prum ruhfel salichos dosreis andersen rayfinned spiralia; do
   python3 scripts/infer_trees.py $NAME $TREE $REP $OUTPATH $ALNPATH processed_model_selection/quantile_model_selection_simulation_control.csv 8
done


