#!/bin/bash


TYPE=$1
python3 selected_models_to_csv.py $TYPE ## -> all_model_selection<pandit/simulation>.csv
Rscript parse_selected_models.R $TYPE  ## -> quantile_model_selection<pandit/simulation>.csv
