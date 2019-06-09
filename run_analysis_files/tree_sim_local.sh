#!/bin/bash
source ~/.bash_profile


TOPPATH=/home/spielman/dms_modelselection/
ALNPATH=${TOPPATH}/simulations/alignments/
OUTPATH=${TOPPATH}/fitted_trees_simulation/
NAME=$1
cd $TOPPATH
for TREE in opisthokonta prum ruhfel salichos dosreis andersen rayfinned spiralia; do
    for REP in {1..20}; do
        python3 run_analysis_files/infer_trees.py $NAME $TREE $REP $OUTPATH $ALNPATH processed_model_selection/quantile_model_selection_simulation.csv 4 
    done
done

