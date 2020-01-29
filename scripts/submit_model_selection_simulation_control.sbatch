#!/bin/bash
#SBATCH --job-name=modelselection
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=100
#SBATCH -p csm
#SBATCH --array=1-20

REP=${SLURM_ARRAY_TASK_ID}
TOPPATH=/csm_data/spielman_lab/aa_phylo_fit_topology
ALNPATH=${TOPPATH}/simulations/alignments/wag_control
OUTPATH=${TOPPATH}/selected_models_simulation_control
mkdir -p $OUTPATH
for NAME in LAC HA NP HIV; do
    for TREE in opisthokonta prum ruhfel salichos dosreis andersen rayfinned spiralia; do
       X=${NAME}_${TREE}_rep${REP}_AA
       if [ ! -f $OUTPATH/${X}.model_selection_log ]; then
         iqtree -nt 8 -m TESTONLY -s ${ALNPATH}/${X}.fasta -st AA
         mv ${ALNPATH}/${X}.fasta.log $OUTPATH/${X}.model_selection_log
         rm $ALNPATH/${X}.fasta.*
        fi
    done
done

