if (!(exists("LOADED"))) source("load.R")
if (exists("LOADED") && LOADED == FALSE) source("load.R") 

############################################################################################
####################################### TABLES #############################################
############################################################################################

#### m1-m5 models for simulations
msel_simulation %>% mutate(simtype = "MutSel") -> mutsel
msel_control %>% 
  mutate(simtype = "control") %>%
  bind_rows(mutsel) %>%
  rowwise() %>%
  mutate(modelm = paste0("m", modelm)) %>%
  count(modelm, name, tree, model_name, simtype) %>%
  mutate(tree  = factor(tree, levels=tree_levels, labels = tree_labels),
         name  = factor(name, levels=name_levels)) %>%
  dplyr::select(simtype, modelm, name, tree, model_name, n) %>%
  arrange(simtype, modelm, name, tree, n) %>%
  rename(sim_replicate = n) %>%
  write_csv(paste0(si_figure_directory, "table_S1.csv"))


##### Simulations where RF = 0 
simulation_rf_fit %>% 
  filter(rf_true_norm == 0) %>%
  count(name_levels, tree_levels, model_levels) %>%
  mutate(simtype = "MutSel") -> sim_rf0
control_rf_fit %>% 
  filter(rf_true_norm == 0) %>%
  count(name_levels, tree_levels, model_levels) %>%
  mutate(simtype = "control") %>%
  bind_rows(sim_rf0) %>%
  rename(name = name_levels, 
         tree = tree_levels, 
         modelm = model_levels,
         number_of_replicates = n) %>%
  write_csv(paste0(si_figure_directory, "table_S2.csv"))

#### Simulations with a SIGNIFICANT AU test
simulation_topology %>%
  filter(whichtest == "au") %>%
  pivot_longer(m1:true, names_to="model", values_to = "pvalue") %>%
  filter(pvalue <= SIG.ALPHA) %>%
  count(model, name, tree) %>%
  arrange(model, name, tree) %>%
  mutate(simtype = "MutSel") -> sim_au
control_topology %>%
  filter(whichtest == "au") %>%
  pivot_longer(m1:true, names_to="model", values_to = "pvalue") %>%
  filter(pvalue <= SIG.ALPHA) %>%
  count(model, name, tree) %>%
  arrange(model, name, tree) %>%
  mutate(simtype = "control") %>%
  bind_rows(sim_au) %>%
  rename(modelm = model,
        sim_replicate = n) %>%
  write_csv(paste0(si_figure_directory, "table_S3.csv"))

#### m1-m5 models for pandit
msel_pandit %>%
  rowwise() %>%
  mutate(modelm = paste0("m", modelm)) %>%
  count(modelm, model_name) %>%
  rename(model = model_name, 
         number_alignments = n) %>%
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
  theme(strip.text = element_text(size=7), 
        legend.key.size = unit(10, "pt"), 
        legend.position = "bottom",
        legend.text = element_text(size=8), 
        legend.title = element_text(size=9), 
        axis.text = element_text(size=7), 
        axis.title = element_text(size=8)) +
  guides(fill = guide_legend(nrow=1)) -> model_fit_simulation_bars_mutsel

control_rf_fit %>%
  ggplot(aes(x = factor(ic.rank), fill = model_levels)) + 
  geom_bar(color="black", size=.2) + 
  facet_grid(name_levels~tree_levels) +
  panel_border() + 
  scale_fill_manual(values=model_colors, name = "Protein Model") +
  xlab("Model rank by BIC") + ylab("Count") + 
  theme(strip.text = element_text(size=7), 
        axis.text = element_text(size=7), 
        axis.title = element_text(size=8)) -> model_fit_simulation_bars_control

l <- get_legend(model_fit_simulation_bars_mutsel)
both <- plot_grid(model_fit_simulation_bars_mutsel + theme(legend.position = "none"),
                  model_fit_simulation_bars_control + theme(legend.position = "none"), 
                  ncol=1, labels = "auto")
both_l <- plot_grid(both, l, ncol=1, rel_heights=c(1,0.05))
save_plot(paste0(si_figure_directory,"model_fit_simulation_bars.pdf"), both_l, base_width=8,base_height=6.25)


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


## NUMBER of false positive nodes

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

ufb_fact_classif %>% 
  ggplot(aes(x = model_levels, y = percent_fp, fill = model_levels)) + 
  geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5, color="grey10")  +     
  scale_fill_manual(values=model_colors, name = "Protein Model") +         
  facet_grid(name_levels~tree_levels, scales="free_y") +
  background_grid() + 
  panel_border() +
  theme(legend.position = "none",
        strip.text = element_text(size=8), 
        axis.text.x = element_text(size=8),
        panel.spacing = unit(0.2, "cm")) +
  xlab("Protein Models") + ylab("Percentage of FP nodes")   -> fp_raw_all


save_plot(paste0(si_figure_directory,"ufb_fpr_all.pdf"), fpr_all, base_width = 12, base_height=10)
save_plot(paste0(si_figure_directory,"ufb_acc_all.pdf"), acc_all, base_width = 12, base_height=10)
save_plot(paste0(si_figure_directory,"ufb_fp_raw_all.pdf"), fp_raw_all, base_width = 12, base_height=10)

  

###################################################################################################
############################## CONTROL SIMULATION FIGURES #########################################
###################################################################################################

################ nRF Boxplots for simulations #################
control_rf_fit %>%
  ggplot(aes(x = model_levels, y = rf_true_norm)) + 
  geom_boxplot(aes(fill = name_levels), outlier.size = 0.3, size=0.25,  width=0.9) + 
  facet_wrap(~tree_levels, scales = "free_y", nrow=2) +
  scale_fill_brewer(palette = "Dark2", name = "Simulation Set") +
  panel_border() +
  background_grid() + 
  xlab("Protein model") + ylab("Normalized Robinson-Foulds Distance") + 
  theme(legend.position = "bottom") -> control_rf_boxplot  
save_plot(paste0(si_figure_directory,"control_rf_boxplot.pdf"), control_rf_boxplot, base_width=10, base_height=4)


control_ufb_fact_classif %>%
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
  geom_hline(yintercept = 0.05, color = "dodgerblue3")  -> fpr_all_control

control_ufb_fact_classif %>%
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
  xlab("Protein Models") + ylab("Accuracy")   -> acc_all_control

control_ufb_fact_classif %>% 
  ggplot(aes(x = model_levels, y = percent_fp, fill = model_levels)) + 
  geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5, color="grey10")  +  
  scale_fill_manual(values=model_colors, name = "Protein Model") +         
  facet_grid(name_levels~tree_levels, scales="free_y") +
  background_grid() + 
  panel_border() +
  theme(legend.position = "none",
        strip.text = element_text(size=8), 
        axis.text.x = element_text(size=8),
        panel.spacing = unit(0.2, "cm")) +
  xlab("Protein Models") + ylab("Percentage of FP nodes")   -> control_fp_raw_all
save_plot(paste0(si_figure_directory,"fp_raw_all.pdf"), fp_raw_all_control, base_width = 12, base_height=10)


save_plot(paste0(si_figure_directory,"control_ufb_fpr_all.pdf"), fpr_all_control, base_width = 12, base_height=10)
save_plot(paste0(si_figure_directory,"control_ufb_acc_all.pdf"), acc_all_control, base_width = 12, base_height=10)
save_plot(paste0(si_figure_directory,"control_ufb_fp_raw_all.pdf"), control_fp_raw_all, base_width = 12, base_height=10)
