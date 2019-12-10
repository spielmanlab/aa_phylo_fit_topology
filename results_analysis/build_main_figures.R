if (is.null(LOADED)) source("load.R") 

#################################################################################
################### Figures associated with simulations #########################
#################################################################################

######### Example BIC distribution to be embedded into flowchart ###########
scheme_name <- "HA"
scheme_tree <- "dosreis"
scheme_rep <- 1
all_sel %>% 
  dplyr::select(name, tree, model, bic, repl) %>%
  filter(name == scheme_name, tree == scheme_tree, repl == scheme_rep) -> data

msel_simulation %>%
  filter(name == scheme_name, tree == scheme_tree, repl == scheme_rep) %>%
  rowwise() %>%
  mutate(modelm2 = paste0("m", modelm),
         model_levels_matrix = paste0(modelm2, " (", model_name, ")")) %>%
  distinct() %>%
  mutate(model_levels= factor(modelm2, levels=c("m1", "m2", "m3", "m4", "m5"))) -> chunks

ggplot(data, aes(x = "", y = bic)) + 
  geom_violin(fill="dodgerblue3", color = "dodgerblue4", alpha=0.3) + 
  geom_point(alpha = 0.1, size=1.5) +
  geom_point(data = chunks, shape=21, aes(x = 1, y = bic, fill = model_levels), size=4) + 
  geom_label(data = chunks,  x = 1.02, hjust="outward", color = "black", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=4, fontface = "bold", alpha=0.8) + 
  scale_fill_manual(values = c(  lighten(lighten(model_colors[1])), model_colors[2:7])) +
  theme(legend.position = "none", axis.ticks.x = element_blank()) +
  xlab("") + ylab("BIC values across all tested models")  -> quant_scheme_plot
ggsave(paste0(figure_directory, "bic_dist_ha_dosreis_ONLYm1-5.pdf"), quant_scheme_plot, width = 3.5, height=4)




############# Barplot showing the m1 matrix for simulations ###############
msel_simulation %>%
  filter(modelm == 1) %>%
  rowwise() %>%
  mutate(model_matrix = str_split(model_name, "\\+")[[1]][1]) %>%
  group_by(name, tree, model_matrix) %>%
  tally() %>% 
  mutate(tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_abbr),
         name_levels  = factor(name, levels=name_levels),
         model_matrix_levels = fct_reorder(model_matrix, n)) %>% 
  ggplot(aes(x = tree_levels, y = n, fill = model_matrix_levels)) + 
  geom_bar(stat="identity") +
  facet_wrap(~name_levels, nrow=2) + scale_fill_brewer(palette = "Dark2", name = paste("m1 Model Matrix")) +
  xlab("Simulation tree") + 
  ylab("Count") + 
  theme(axis.text.x = element_text(size=8)) -> selected_m1_plot    
save_plot(paste0(figure_directory, "selected_m1_simulation_barplot.pdf"), selected_m1_plot, base_width=10, base_height = 4) 


################ nRF Boxplots for simulations #################
simulation_rf_fit %>%
  ggplot(aes(x = model_levels, y = rf_true_norm)) + 
  geom_boxplot(aes(fill = name_levels), outlier.size = 0.3, size=0.25,  width=0.9) + 
  facet_wrap(~tree_levels, scales = "free_y", nrow=2) +
  scale_fill_brewer(palette = "Set2", name = "Simulation Set") +
  panel_border() +
  background_grid() + 
  xlab("Protein model") + ylab("Normalized Robinson-Foulds Distance") + 
  theme(legend.position = "bottom") -> simulation_rf_boxplot  
save_plot(paste0(figure_directory,"simulation_rf_boxplot.pdf"), simulation_rf_boxplot, base_width=10, base_height=4)



######################## FPR and accuracy for HA as example ####################
demo_name <- "HA"
ufb_fact_classif %>%
  filter(name_levels == demo_name) %>%
  ggplot(aes(x = model_levels, y = FPR, fill = model_levels)) + 
  geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5, color="grey10")  +     
  scale_fill_manual(values=model_colors, name = "Protein Model") +         
  facet_grid(~tree_levels, scales="free_y") +
  background_grid() + 
  panel_border() +
  theme(legend.position = "none",
        strip.text = element_text(size=8), 
        axis.text.x = element_text(size=7),
        panel.spacing = unit(0.2, "cm")) +
  xlab("Protein Models") + ylab("False positive rate") +
  geom_hline(yintercept = 0.05, color = "dodgerblue3")  -> fpr_ha

ufb_fact_classif %>%
  filter(name_levels == demo_name) %>%
  ggplot(aes(x = model_levels, y = accuracy, fill = model_levels)) + 
  geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5, color="grey10")  +     
  scale_fill_manual(values=model_colors, name = "Protein Model") +         
  facet_grid(~tree_levels, scales="free_y") +
  background_grid() + 
  panel_border() +
  theme(legend.position = "none",
        strip.text = element_text(size=8), 
        axis.text.x = element_text(size=7),
        panel.spacing = unit(0.2, "cm")) +
  xlab("Protein Models") + ylab("Accuracy")   -> acc_ha

ufb_ha_grid <- plot_grid(fpr_ha, acc_ha, nrow=2, labels= "auto", scale=0.98)
save_plot(paste0(figure_directory,"ufb_ha_grid.pdf"), ufb_ha_grid, base_width = 12, base_height=4)


#################################################################################
################### Figures associated with PANDIT analysis #####################
#################################################################################


############# Pandit model fit: a) m1 models, b) JC/GTR in BIC ranks, c) when does GTR *not* overfit? #########
msel_pandit %>%
  filter(modelm == 1) %>%
  rowwise() %>%
  mutate(model_matrix = str_split(model_name, "\\+")[[1]][1]) %>%
  group_by(model_matrix) %>%
  tally() %>% 
  mutate(model_matrix_levels = fct_reorder(model_matrix, n)) %>%
  ggplot(aes(x = model_matrix_levels, y = n), fill = "grey60") + 
  geom_col() +
  geom_text(aes(x = model_matrix_levels, y = n+4, label = n), size=3)+
  xlab("m1 Model Matrix") + 
  ylab("Count") +
  theme(axis.text.x = element_text(size=7))-> m1_pandit_model_plot



pandit_ranks %>%    
  filter(model %in% c("GTR20", "poisson")) %>%
  count(model, ic.rank) %>%
  mutate(model_levels = factor(model, levels=c("GTR20", "poisson"), labels = c("GTR", "JC"))) %>%
  ggplot(aes(x = factor(ic.rank), fill = model_levels, y = n)) + 
  geom_col(color="black", size=.2, position = position_dodge2(width=1, preserve="single")) + 
  geom_text(aes(x = factor(ic.rank), y = n+6, label = n), size=2.25, position = position_dodge(width=1)) + 
  scale_fill_brewer(palette = "Set2", name = "") +
  xlab("Overall model rank") + ylab("Count") -> gtr_jc_pandit_bars




pandit_ranks %>% 
  filter(model == "GTR20") %>% 
  ggplot(aes(x = ic.rank, y = ntaxa, color = tl)) + 
  geom_point(alpha=0.7, size=2) + geom_smooth(method = "lm", color= "grey20") + 
  scale_x_continuous(breaks=1:6) + 
  scale_color_continuous(name = "Treelength") + 
  xlab("GTR model rank") + ylab("Number of taxa") + 
  annotate("text", x = 5, y=325, label = "R^2 == 0.61", parse=TRUE) -> gtr20_rank_plot

plot_grid(m1_pandit_model_plot, gtr_jc_pandit_bars, gtr20_rank_plot, scale=c(0.93, 0.97, 0.93), nrow=1, labels="auto") -> pandit_model_grid
save_plot(paste0(figure_directory, "pandit_model_bars.pdf"), pandit_model_grid, base_width=12, base_height=2.5)



################### Ridgeplot of PANDIT nRF all-to-all comparison ##################
pandit_rf %>%
  group_by(model1, model2) %>% 
  summarize(medianrf = median(rf)) %>% ## 0.294 is min
  left_join(pandit_rf) %>%
  ungroup() %>%
  mutate(model1 = factor(model1, levels=model_levels, labels = model_labels), 
         model2 = factor(model2, levels=model_levels, labels = model_labels)) %>%
  mutate(hasm1 = (model1 == "m1" | model2 == "m1")) %>%
  unite(model_pair, model1, model2, sep = " - ") %>%
  mutate(model_pair_fct = fct_reorder(model_pair, rf, .desc=TRUE)) -> pandit_rf_ridgedata

pandit_rf_ridgedata %>% 
  dplyr::select(model_pair_fct, hasm1) %>% 
  distinct() %>% 
  arrange(model_pair_fct) %>% 
  mutate(face = ifelse(hasm1, "bold", "plain")) %>% 
  pull(face) -> font_face
pandit_rf_ridgedata %>% 
  ggplot(aes(x = rf, y = fct_reorder(model_pair, rf, .desc=TRUE), fill = medianrf)) + 
  geom_density_ridges(scale=1.4, quantile_lines=TRUE, quantiles=2, rel_min_height = 0.01) +
  ylab("") + xlab("Normalized Robinson-Foulds Distance") + 
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), position="right") +
  background_grid(colour.major = "grey70")+
  scale_fill_gradient(low = "#FFF68F", high = "dodgerblue3", name = "Median nRF",
                      guide = guide_colorbar(
                        direction = "horizontal",
                        label.position = "bottom",
                        title.position = "left",
                        barwidth = grid::unit(1.5, "in"),
                        barheight = grid::unit(0.2, "in")
                      )
  ) +
  theme(
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_text(size=11, face = font_face),
    axis.text.x = element_text(size=11)
  )  -> pandit_rf_plot


save_plot(paste0(figure_directory, "pandit_ridgeplot.pdf"), pandit_rf_plot, base_width=5, base_height=7)


################### PANDIT AU test results ##################
pandit_topology %>%
  filter(whichtest == "au") %>%
  mutate(totaln = n()) %>%
  gather(model, pvalue, m1:GTR20) %>% 
  mutate(notsig = pvalue >= SIG.ALPHA) %>%
  count(model, notsig, totaln) %>%
  mutate(perc_in_conf=n/totaln) %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) -> pandit_au_summary

pandit_au_summary  %>% filter(notsig == TRUE) %>% mutate(label = paste0(100*(perc_in_conf), "%")) -> pandit_au_summary_false
pandit_au_summary %>%
  mutate(notsig = factor(notsig, levels=c("TRUE", "FALSE"))) %>%
  ggplot(aes(x = model_levels, y = perc_in_conf, fill = notsig)) + 
  geom_col() +
  geom_text(data = pandit_au_summary_false, aes(x = model_levels, y = 1.04, label = label), size=3) + 
  scale_fill_brewer(palette = "Set2", name = "In m1 confidence set") +
  xlab("Protein models") + ylab("Proportion of alignments") +
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(nrow=1, title.position = "left")) -> pandit_au_barplot


pandit_topology %>% filter(whichtest == "au", m4 >=SIG.ALPHA, m5>=SIG.ALPHA, poisson <SIG.ALPHA) %>% nrow()  -> poisson_only   
pandit_topology %>% filter(whichtest == "au", m4 >=SIG.ALPHA, m5<SIG.ALPHA, poisson >=SIG.ALPHA) %>% nrow()  -> m5_only
pandit_topology %>% filter(whichtest == "au", m4 <SIG.ALPHA,  m5>=SIG.ALPHA, poisson >=SIG.ALPHA) %>% nrow() -> m4_only
pandit_topology %>% filter(whichtest == "au", m4 >=SIG.ALPHA, m5<SIG.ALPHA, poisson <SIG.ALPHA) %>% nrow()   -> m5_poisson       
pandit_topology %>% filter(whichtest == "au", m4 <SIG.ALPHA,  m5<SIG.ALPHA, poisson >=SIG.ALPHA) %>% nrow()  -> m4_m5
pandit_topology %>% filter(whichtest == "au", m4 <SIG.ALPHA,  m5>=SIG.ALPHA, poisson <SIG.ALPHA) %>% nrow()  -> m4_poisson
pandit_topology %>% filter(whichtest == "au", m4 <SIG.ALPHA,  m5<SIG.ALPHA, poisson <SIG.ALPHA) %>% nrow()   -> all_three        
venn_numbers <- tibble(type = c("m4_only", "m5_only", "poisson_only", "m4_m5", "m4_poisson", "m5_poisson", "all"), 
                       labels = c(m4_only, m5_only, poisson_only, m4_m5, m4_poisson, m5_poisson, all_three),                       x = c(0,     1.4,  -1.4,    0.75,   -0.75,   0,   0),   
                       y = c(1.5,  -0.7,   -0.7,    0.5,    0.5,    -1,   0))                                
df.venn <- data.frame(x = c(0, 0.866, -0.866),
                      y = c(1, -0.5, -0.5),
                      labels = c('m4', 'm5', 'JC'))
ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = 0.6, size = 1, colour = 'grey50') +
  theme_void() + 
  coord_fixed(clip = "off") + 
  theme(legend.position = 'none') +
  scale_fill_brewer(palette = "Set2") +
  labs(fill = NULL) +
  annotate("text", x = venn_numbers$x, y = venn_numbers$y, label = venn_numbers$labels, size = 5) +
  annotate("text", x = 0, y = 2.65, label = "m4",fontface = "bold", size=4) +
  annotate("text", x = -2.3, y = -1.5,  label = "JC",fontface = "bold", size=4) +
  annotate("text", x = 2.3, y = -1.5, label = "m5",fontface = "bold", size=4) -> disagree_m1_venn

pandit_au_plot <- plot_grid(pandit_au_barplot, disagree_m1_venn, labels="auto", scale=c(0.95, 0.8))
save_plot(paste0(figure_directory,"pandit_au.pdf"), pandit_au_plot, base_width=8, base_height=3.5)


