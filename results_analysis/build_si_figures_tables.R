if (!(exists("LOADED"))) source("load.R")
if (exists("LOADED") && LOADED == FALSE) source("load.R") 

############################################################################################
####################################### TABLES #############################################
############################################################################################

#### m1-m5 models for simulations
msel_simulation %>%
  rowwise() %>%
  mutate(modelm = paste0("m", modelm)) %>%
  count(modelm, name, tree, model_name) %>%
  mutate(tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels),
         name_levels  = factor(name, levels=name_levels)) %>%
  dplyr::select(modelm, name_levels, tree_levels, model_name, n) %>%
  rename("Model rank" = modelm, 
         "Model name" = model_name,
         "Simulation DMS" = name_levels,
         "Simulation Tree" = tree_levels,
         "Number of replicates" = n) %>%
  write_csv(paste0(si_figure_directory, "table_S1.csv"))


##### Simulations where RF = 0 
simulation_rf_fit %>% 
  filter(rf_true_norm == 0) %>%
  count(name_levels, tree_levels, model_levels) %>% 
  write_csv(paste0(si_figure_directory, "table_S2.csv"))

#### Simulations with a SIGNIFICANT AU test
simulation_topology %>%
  filter(whichtest == "au") %>%
  pivot_longer(m1:true, names_to="model", values_to = "pvalue") %>%
  filter(pvalue <= SIG.ALPHA) %>%
  count(model, name, tree) %>%
  arrange(model, name, tree) %>%
  write_csv(paste0(si_figure_directory, "table_S3.csv"))

#### m1-m5 models for pandit
msel_pandit %>%
  rowwise() %>%
  mutate(modelm = paste0("m", modelm)) %>%
  count(modelm, model_name) %>%
  rename("Model rank" = modelm, 
         "Model name" = model_name, 
         "Number of alignments" = n) %>%
  write_csv(paste0(si_figure_directory, "table_S4.csv"))


############################################################################################
####################################### FIGURES ############################################
############################################################################################


########## Barplots for where JC and GTR fall in the scheme of model selection #########
simulation_rf_fit %>%
  ggplot(aes(x = factor(ic.rank), fill = model_levels)) + 
  geom_bar(color="black", size=.2) + 
  facet_grid(name_levels~tree_levels) +
  panel_border() + 
  scale_fill_manual(values=model_colors, name = "Protein Model") +
  xlab("Model rank by BIC") + ylab("Count") + 
  theme(strip.text = element_text(size=8), 
        legend.position = "bottom", 
        legend.key.size = unit(6, "pt"), 
        legend.text = element_text(size=7), 
        legend.title = element_text(size=7), 
        axis.text = element_text(size=7), 
        axis.title = element_text(size=8)) +
  guides(fill = guide_legend(nrow=1)) -> model_fit_simulation_bars
l <- get_legend(model_fit_simulation_bars)
both <- plot_grid(model_fit_simulation_bars + theme(legend.position = "none"), l, nrow=2, rel_heights=c(1,0.05))
save_plot(paste0(si_figure_directory,"model_fit_simulation_bars.pdf"), both, base_width=10)


################ Protein model correlations #################
m1_models <- model_comp %>% filter(model1 %in% c("JTT", "HIVb", "WAG"), model2 %in% c("JTT", "HIVb", "WAG"))
ggplot(model_comp, aes(x = fct_relevel(model1, c("JTT", "HIVb", "WAG")), y = fct_relevel(model2, c("JTT", "HIVb", "WAG")), fill = r)) + 
  geom_tile() + 
  xlab("") + ylab("") + 
  scale_fill_gradient(name = "Pearson Correlation Coefficient", low = "red", high = "yellow") +
  theme(axis.text.x = element_text(size=8), legend.position = "bottom") + 
  geom_tile(data = m1_models, aes(x = model1, y = model2), color = "black", size=0.8) +
  geom_text(aes(label = round(r, 3)), size=2) -> all_models_r
save_plot(paste0(si_figure_directory, "model_pearson_heatmap.pdf"), all_models_r, base_width=10, base_height=6)


##################### Entropy of the simulations ##########################
ggplot(entropy, aes(x = name_levels, y = total_entropy, color = name_levels)) + 
  geom_sina(size=0.5) + 
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~tree_levels, nrow=2) + 
  stat_summary(geom = "point", color = "black", size=2, pch=18) + 
  xlab("Simulations") + ylab("Alignment entropy") + 
  panel_border() + theme(legend.position = "none") -> entropy_sina
save_plot(paste0(si_figure_directory, "entropy_simulations_sina.pdf"), entropy_sina, base_width=8, base_height=3)



#################### FPR and Accuracy for ALL simulations ##################
ufb_fact_classif %>%
  ggplot(aes(x = model_levels, y = FPR, fill = model_levels)) + 
  geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5, color="grey10")  +     
  scale_fill_manual(values=model_colors, name = "Protein Model") +         
  facet_grid(name_levels~tree_levels, scales="free_y") +
  background_grid() + 
  panel_border() +
  theme(legend.position = "none",
        strip.text = element_text(size=8), 
        axis.text.x = element_text(size=8),
        panel.spacing = unit(0.2, "cm")) +
  xlab("Protein Models") + ylab("False positive rate") +
  geom_hline(yintercept = 0.05, color = "dodgerblue3")  -> fpr_all

ufb_fact_classif %>%
  ggplot(aes(x = model_levels, y = accuracy, fill = model_levels)) + 
  geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5, color="grey10")  +     
  scale_fill_manual(values=model_colors, name = "Protein Model") +         
  facet_grid(name_levels~tree_levels, scales="free_y") +
  background_grid() + 
  panel_border() +
  theme(legend.position = "none",
        strip.text = element_text(size=8), 
        axis.text.x = element_text(size=8),
        panel.spacing = unit(0.2, "cm")) +
  xlab("Protein Models") + ylab("Accuracy")   -> acc_all
save_plot(paste0(si_figure_directory,"ufb_fpr_all.pdf"), fpr_all, base_width = 12, base_height=10)
save_plot(paste0(si_figure_directory,"ufb_acc_all.pdf"), acc_all, base_width = 12, base_height=10)

###################################################################################################
########################### MISCELLANEOUS EYEBALLING OF SOME DATA #################################
###################################################################################################


### which specific simulations had RF=0?
# simulation_rf_fit %>% 
#   filter(rf_true_norm == 0) %>%
#   count(name_levels, tree_levels, model_levels) %>%
#   print.data.frame()
# name_levels  tree_levels model_levels  n
# 1           NP     Spiralia          GTR  1
# 2           NP Opisthokonta           m1  1
# 3           NP Opisthokonta           m3  1
# 4           NP Opisthokonta           JC  1
# 5           NP Opisthokonta          GTR  3
# 6           NP        Yeast           m1  4
# 7           NP        Yeast           m2  3
# 8           NP        Yeast           m3  1
# 9           NP        Yeast           m4  3
# 10          NP        Yeast           JC  1
# 11          NP        Yeast          GTR  4
# 12          HA     Spiralia           m2  1
# 13          HA     Spiralia           m3  1
# 14          HA     Spiralia          GTR  1
# 15          HA Opisthokonta           m2  1
# 16          HA Opisthokonta           m3  1
# 17          HA Opisthokonta          GTR  4
# 18          HA        Yeast           m1  2
# 19          HA        Yeast           m2  2
# 20          HA        Yeast           m4  1
# 21          HA        Yeast           m5  1
# 22          HA        Yeast           JC  1
# 23          HA        Yeast          GTR  4
# 24         HIV     Spiralia           m1  3
# 25         HIV     Spiralia           m2  2
# 26         HIV     Spiralia           m3  3
# 27         HIV     Spiralia           m4  4
# 28         HIV     Spiralia           JC  4
# 29         HIV     Spiralia          GTR  5
# 30         HIV Opisthokonta           m1  9
# 31         HIV Opisthokonta           m2  9
# 32         HIV Opisthokonta           m3  4
# 33         HIV Opisthokonta           m4  4
# 34         HIV Opisthokonta           JC  6
# 35         HIV Opisthokonta          GTR 16
# 36         HIV        Yeast           m1  2
# 37         HIV        Yeast           m2  2
# 38         HIV        Yeast           m3  2
# 39         HIV        Yeast           m4  2
# 40         HIV        Yeast           m5  1
# 41         HIV        Yeast           JC  3
# 42         HIV        Yeast          GTR  6

### which trees overall showed SIGNIFICANT AU tests?
# simulation_topology %>%
#   filter(whichtest == "au") %>%
#   pivot_longer(m1:true, names_to="model", values_to = "pvalue") %>%
#   filter(pvalue <= SIG.ALPHA) %>%
#   count(model)
# A tibble: 3 x 2
# model     n
# <chr> <int>
# 1 GTR20     1
# 2 m5       25
# 3 true     72

### which specific simulations showed SIGNIFICANT AU tests?
# simulation_topology %>%
#   filter(whichtest == "au") %>%
#   pivot_longer(m1:true, names_to="model", values_to = "pvalue") %>%
#   filter(pvalue <= SIG.ALPHA) %>%
#   count(model, name, tree) %>% 
#   arrange(model, name, tree) %>%
#   print.data.frame()
# model name         tree  n
# 1  GTR20  HIV      dosreis  1
# 2     m5   HA     andersen  1
# 3     m5   HA opisthokonta  1
# 4     m5  HIV opisthokonta  5
# 5     m5  HIV       ruhfel  4
# 6     m5  HIV     salichos  1
# 7     m5  HIV     spiralia  1
# 8     m5  LAC      dosreis  1
# 9     m5  LAC opisthokonta  3
# 10    m5  LAC       ruhfel  1
# 11    m5  LAC     salichos  1
# 12    m5   NP      dosreis  1
# 13    m5   NP opisthokonta  2
# 14    m5   NP       ruhfel  3
# 15  true   HA     andersen  5
# 16  true   HA opisthokonta  1
# 17  true   HA         prum  9
# 18  true   HA    rayfinned  2
# 19  true   HA       ruhfel  3
# 20  true  HIV         prum  1
# 21  true  LAC     andersen 13
# 22  true  LAC      dosreis  2
# 23  true  LAC         prum  8
# 24  true  LAC    rayfinned 11
# 25  true  LAC       ruhfel  6
# 26  true   NP     andersen  4
# 27  true   NP         prum  6
# 28  true   NP       ruhfel  1

