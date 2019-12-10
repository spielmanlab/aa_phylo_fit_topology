This directory contains scripts and (most) CSV files for analyzing, modeling, visualizing results.

Contents:

+ `csv_files/` contains most summary files for reading, processing. Produced by scripts in `../scripts/`
+ `linear_model_analysis.Rmd` does the linear models presented in main text.
+ `build_main_figures.R` builds all the figures in the main manuscript, which get saved to `main_figures/`
+ `build_si_tables_figures.R` builds all the figures and tables for SI, which get saved to `si_figures_tables/`
+ `load.R` is called by all R scripts (including Rmd) for loading up libraries, loading and formatting data
