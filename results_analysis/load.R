################################################################
#                                                              #    
# Load all libraries, data for generating tables and figures.  #
#                                                              #
################################################################
LOADED <- TRUE

library(tidyverse)
library(cowplot)
library(xtable)
library(colorspace)
library(ggridges)
library(ggforce)

theme_set(theme_classic() + theme(axis.line = element_line(colour = "grey10"),
                                  strip.background = element_rect(size=0.5, fill="grey90")))
data_path <- "csv_files/"
figure_directory <- "main_figures/"
si_figure_directory <- "si_figures_tables/"
SIG.ALPHA <- 0.01 ## significance threshold
########################## Factor levels and labeling ################################ 
name_levels <- c("LAC", "NP", "HA", "HIV")

model_levels <- c("m1", "m2", "m3", "m4", "m5", "poisson", "GTR20")
model_labels <- c("m1", "m2", "m3", "m4", "m5", "JC", "GTR")
m1to5_cols <- sequential_hcl(5, palette = "ylorrd")
poisson_col <- "grey40"
gtr20_col   <- "grey80" 
model_colors <- c(m1to5_cols, poisson_col, gtr20_col)
model_colors_nom1 <- c(m1to5_cols[2:5], poisson_col, gtr20_col)

model_levels_nom1 <- c("m2", "m3", "m4", "m5", "poisson", "GTR20")
model_labels_nom1 <- c("m2", "m3", "m4",  "m5", "JC", "GTR")

tree_levels <- c("ruhfel", "rayfinned", "dosreis", "prum", "andersen", "spiralia", "opisthokonta", "salichos")
tree_labels <- c("Green Plant", "Ray-finned fish", "Mammals", "Aves", "Lassa Virus", "Spiralia", "Opisthokonta", "Yeast")
tree_labels_abbr <- c("Plant", "Fish", "Mammals", "Aves", "Lassa", "Spiralia", "Opis.", "Yeast")

################################### Read in all data ##################################

### basic information about datasets
sim_info    <- read_csv("../simulations/simulation_trees_ntaxa.csv") 
pandit_info <- read_csv("../pandit_aa_alignments/pandit_info.csv")

### RF and BIC results
simulation_rf_fit  <- read_csv(paste0(data_path, "rf_fit_simulation.csv"),guess_max = 10000)
simulation_rf_fit %<>% 
    left_join(sim_info) %>%
    mutate(max_rf = 2*ntaxa - 6) %>%  ### twice the number of internal edges
    mutate(rf_true_norm = rf/max_rf) %>%
    dplyr::select(-AIC, -AICc, -k,-logl, -rf, -treelength, -max_rf, -ntaxa) %>%
    distinct() %>%
    group_by(name, tree, rep) %>%
    mutate(ic.rank = as.integer(rank(BIC))) %>%
    ungroup() %>%
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
           tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels),
           name_levels  = factor(name, levels=name_levels)) %>%
    dplyr::select(-model, -tree, -name)

pandit_rf  <- read_csv(paste0(data_path, "rf_pandit.csv"), guess_max = 10000) 
pandit_fit <- read_csv(paste0(data_path, "fit_pandit.csv"), guess_max = 10000) 
pandit_rf %>% 
    left_join(pandit_info) %>%
    mutate(max_rf = 2*ntaxa - 6) %>%
    mutate(rf     = rf/max_rf) -> pandit_rf
pandit_fit %>% left_join(pandit_info) -> pandit_fit
pandit_fit %>%
    group_by(name) %>%
    mutate(ic.rank = as.integer(rank(BIC))) %>%
    left_join(pandit_info) %>%
    ungroup() -> pandit_ranks

### Topology test results
simulation_topology <- read_csv(paste0(data_path, "topology_tests_simulation.csv"))
pandit_topology     <- read_csv(paste0(data_path, "topology_tests_pandit.csv"))

## actual selected models
msel_simulation <- read_csv("../processed_model_selection/quantile_model_selection_simulation.csv")
all_sel         <- read_csv("../processed_model_selection/all_model_selection_simulation.csv")
msel_pandit     <- read_csv("../processed_model_selection/quantile_model_selection_pandit.csv")


## Bootstrap analysis
sim_ufb <- read_csv(paste0(data_path, "ufb_splits_simulation.csv"))

sim_ufb %>% 
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
        tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels),
        name_levels  = factor(name, levels=name_levels), 
        rep = factor(rep)) -> ufb_fact
ufb_fact %>% 
    count(model_levels, tree_levels, name_levels, rep, classif) %>%
    pivot_wider(names_from = classif, values_from = n) %>%
    ungroup() %>%
    replace_na(list(FP = 0, FN = 0, TP = 0, TN = 0)) %>%
    mutate(FPR = ifelse(is.nan(FP / (TN+FP)), 0, FP / (TN+FP)), 
           accuracy = (TP + TN)/(TP+TN+FP+FN)) -> ufb_fact_classif


##### entropy (for distinguishing simulations) ######
entropy <- read_csv(paste0(data_path, "simulations_entropy.csv")) %>%
                mutate(name_levels = factor(name, levels=name_levels),
                       tree_levels = factor(tree, levels = tree_levels, labels = tree_labels))

#### diffs among models ####
model_comp <- read_csv(paste0(data_path, "all_models_pearson.csv"))
all_models <- unique(c(model_comp$model1, model_comp$model2))
all_models <- all_models[!all_models %in% c("JTT", "HIVb", "WAG")]
all_models <- c("JTT", "HIVb", "WAG", all_models)
model_comp$model1 <- factor(model_comp$model1, levels = all_models)
model_comp$model2 <- factor(model_comp$model2, levels = all_models)
