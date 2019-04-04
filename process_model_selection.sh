#!/bin/bash

python3 selected_models_to_csv.py ## -> all_model_selection.csv
Rscript parse_selected_models.R   ## -> quantile_model_selection.csv