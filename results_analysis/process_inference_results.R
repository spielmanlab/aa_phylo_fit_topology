source("load.R") # Load all libraries and data
#stop()


############################### Simulation figures ################################

#### simulation fit results
simulation_rf_fit %>%
  ggplot(aes(x = factor(ic.rank), fill = model_levels)) + 
    geom_bar(color="black", size=.2) + 
    facet_grid(name_levels~tree_levels) +
    panel_border() + 
    #scale_fill_brewer(palette = "RdYlBu", name = "Protein Model") +
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
save_plot(paste0(figure_directory,"model_fit_simulation_bars.pdf"), both, base_width=10)
  
  


######## simulation RF boxplots
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


simulation_rf_fit %>% 
    filter(rf_true_norm == 0) %>%
    count(name_levels, tree_levels, model_levels) %>% 
    xtable()
# % latex table generated in R 3.6.1 by xtable 1.8-4 package
# % Thu Dec  5 10:55:07 2019
# \begin{table}[ht]
# \centering
# \begin{tabular}{rlllr}
# \hline
# & name\_levels & tree\_levels & model\_levels & n \\ 
# \hline
# 1 & NP & Spiralia & GTR &   1 \\ 
# 2 & NP & Opisthokonta & m1 &   1 \\ 
# 3 & NP & Opisthokonta & m3 &   1 \\ 
# 4 & NP & Opisthokonta & JC &   1 \\ 
# 5 & NP & Opisthokonta & GTR &   3 \\ 
# 6 & NP & Yeast & m1 &   4 \\ 
# 7 & NP & Yeast & m2 &   3 \\ 
# 8 & NP & Yeast & m3 &   1 \\ 
# 9 & NP & Yeast & m4 &   3 \\ 
# 10 & NP & Yeast & JC &   1 \\ 
# 11 & NP & Yeast & GTR &   4 \\ 
# 12 & HA & Spiralia & m2 &   1 \\ 
# 13 & HA & Spiralia & m3 &   1 \\ 
# 14 & HA & Spiralia & GTR &   1 \\ 
# 15 & HA & Opisthokonta & m2 &   1 \\ 
# 16 & HA & Opisthokonta & m3 &   1 \\ 
# 17 & HA & Opisthokonta & GTR &   4 \\ 
# 18 & HA & Yeast & m1 &   2 \\ 
# 19 & HA & Yeast & m2 &   2 \\ 
# 20 & HA & Yeast & m4 &   1 \\ 
# 21 & HA & Yeast & m5 &   1 \\ 
# 22 & HA & Yeast & JC &   1 \\ 
# 23 & HA & Yeast & GTR &   4 \\ 
# 24 & HIV & Spiralia & m1 &   3 \\ 
# 25 & HIV & Spiralia & m2 &   2 \\ 
# 26 & HIV & Spiralia & m3 &   3 \\ 
# 27 & HIV & Spiralia & m4 &   4 \\ 
# 28 & HIV & Spiralia & JC &   4 \\ 
# 29 & HIV & Spiralia & GTR &   5 \\ 
# 30 & HIV & Opisthokonta & m1 &   9 \\ 
# 31 & HIV & Opisthokonta & m2 &   9 \\ 
# 32 & HIV & Opisthokonta & m3 &   4 \\ 
# 33 & HIV & Opisthokonta & m4 &   4 \\ 
# 34 & HIV & Opisthokonta & JC &   6 \\ 
# 35 & HIV & Opisthokonta & GTR &  16 \\ 
# 36 & HIV & Yeast & m1 &   2 \\ 
# 37 & HIV & Yeast & m2 &   2 \\ 
# 38 & HIV & Yeast & m3 &   2 \\ 
# 39 & HIV & Yeast & m4 &   2 \\ 
# 40 & HIV & Yeast & m5 &   1 \\ 
# 41 & HIV & Yeast & JC &   3 \\ 
# 42 & HIV & Yeast & GTR &   6 \\ 
# \hline
# \end{tabular}
# \end{table}


simulation_rf_fit %>% 
  filter(rf_true_norm == 0) %>%
  count(model_levels) %>%
  arrange(desc(n))
# # A tibble: 7 x 2
# model_levels     n
# <fct>        <int>
# 1 GTR             44
# 2 m1              21
# 3 m2              20
# 4 JC              16
# 5 m4              14
# 6 m3              13
# 7 m5               2

simulation_topology %>%
  filter(whichtest == "au") %>%
  pivot_longer(m1:true, names_to="model", values_to = "pvalue") %>%
  filter(pvalue <= SIG.ALPHA) %>%
  count(model)
# A tibble: 3 x 2
# model     n
# <chr> <int>
# 1 GTR20     1
# 2 m5       25
# 3 true     72

simulation_topology %>%
  filter(whichtest == "au") %>%
  pivot_longer(m1:true, names_to="model", values_to = "pvalue") %>%
  filter(pvalue <= SIG.ALPHA) %>%
  count(model, name, tree) %>% 
  arrange(model, name, tree) %>%
  print.data.frame()
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


demo_name <- "HA"
ufb_fact_classif %>%
    filter(name_levels == demo_name) %>%
    ggplot(aes(x = model_levels, y = FPR, fill = model_levels)) + 
        geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5, color="grey10")  +     
        scale_fill_manual(values=model_colors, name = "Protein Model") +         
        facet_grid(~tree_levels, scales="free_y") +
        #scale_y_continuous(breaks=c(0, 0.02, 0.04, 0.06)) + 
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
        #scale_y_continuous(breaks=c(0, 0.02, 0.04, 0.06)) + 
        background_grid() + 
        panel_border() +
        theme(legend.position = "none",
              strip.text = element_text(size=8), 
              axis.text.x = element_text(size=7),
              panel.spacing = unit(0.2, "cm")) +
        xlab("Protein Models") + ylab("Accuracy")   -> acc_ha


ufb_ha_grid <- plot_grid(fpr_ha, acc_ha, nrow=2, labels= "auto", scale=0.98)
save_plot(paste0(figure_directory,"ufb_ha_grid.pdf"), ufb_ha_grid, base_width = 12, base_height=4)
             
             
             
             
ufb_fact_classif %>%
    ggplot(aes(x = model_levels, y = FPR, fill = model_levels)) + 
        geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5, color="grey10")  +     
        scale_fill_manual(values=model_colors, name = "Protein Model") +         
        facet_grid(name_levels~tree_levels, scales="free_y") +
        #scale_y_continuous(breaks=c(0, 0.02, 0.04, 0.06)) + 
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
        #scale_y_continuous(breaks=c(0, 0.02, 0.04, 0.06)) + 
        background_grid() + 
        panel_border() +
        theme(legend.position = "none",
              strip.text = element_text(size=8), 
              axis.text.x = element_text(size=8),
              panel.spacing = unit(0.2, "cm")) +
        xlab("Protein Models") + ylab("Accuracy")   -> acc_all
save_plot(paste0(figure_directory,"ufb_fpr_all.pdf"), fpr_all, base_width = 12, base_height=10)
save_plot(paste0(figure_directory,"ufb_acc_all.pdf"), acc_all, base_width = 12, base_height=10)


################################## PANDIT figures ###########################################


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

 #### PANDIT fit results
pandit_fit %>%
    group_by(name) %>%
    mutate(ic.rank = as.integer(rank(BIC))) %>%
    left_join(pandit_info) -> pandit_ranks
    
pandit_ranks %>%    
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) %>%
    ggplot(aes(x = factor(ic.rank), fill = model_levels)) + 
    geom_bar(color="black", size=.2) + 
    panel_border() + 
    scale_fill_manual(values=model_colors, name = "") +
    xlab("Model rank by BIC") + ylab("Count")-> model_fit_pandit_bars
    #theme(legend.key.size = unit(6, "pt"), 
    #      legend.text = element_text(size=7), 
    #      legend.title = element_text(size=7))
    #guides( fill = guide_legend(nrow=1))
     
  
pandit_ranks %>% 
    filter(model == "GTR20") %>% 
    ggplot(aes(x = ic.rank, y = ntaxa, color = tl)) + 
        geom_point(alpha=0.7, size=2) + geom_smooth(method = "lm", color= "grey20") + 
        scale_x_continuous(breaks=1:6) + 
        scale_color_continuous(name = "Treelength") + 
        xlab("GTR model rank") + ylab("Number of taxa") + 
        annotate("text", x = 5, y=325, label = "R^2 == 0.61", parse=TRUE) -> gtr20_rank_plot
        
        
        
plot_grid(m1_pandit_model_plot, model_fit_pandit_bars, gtr20_rank_plot, scale=0.93, nrow=1, labels="auto") -> pandit_model_grid

save_plot(paste0(figure_directory, "pandit_model_bars.pdf"), pandit_model_grid, base_width=12, base_height=2.5)
  





pandit_rf %>%
    group_by(model1, model2) %>% 
    summarize(medianrf = median(rf)) %>% ## 0.294is min
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
        )  -> rf_pandit_plot
 
save_plot(paste0(figure_directory, "pandit_ridgeplot.pdf"), rf_pandit_plot, base_width=5, base_height=7)
        
        
    
    

pandit_topology %>%
    filter(whichtest == "au") %>%
    gather(model, pvalue, m1:GTR20) %>% 
    mutate(notsig = pvalue >= SIG.ALPHA) %>%
    group_by(model, notsig) %>%
    tally() %>%
    mutate(perc_in_conf=n/200) %>%
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) -> pandit_au_summary

pandit_au_summary  %>% filter(notsig == TRUE) %>% mutate(label = paste0(100*(perc_in_conf), "%")) -> pandit_au_summary_false
pandit_au_summary %>%
  mutate(notsig = factor(notsig, levels=c("TRUE", "FALSE"))) %>%
    ggplot(aes(x = model_levels, y = perc_in_conf, fill = notsig)) + 
        geom_col() +
        geom_text(data = pandit_au_summary_false, aes(x = model_levels, y = 1.04, label = label), size=3) + 
        #scale_fill_manual(values=c("seagreen3", "coral2"), name = "In M1 confidence set") +
        scale_fill_hue(l=50, name = "In m1 Confidence Set") +
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
  geom_circle(alpha = 0.4, size = 1, colour = 'grey50') +
  theme_void() + 
  coord_fixed(clip = "off") + 
theme(legend.position = 'none') +
  scale_fill_brewer(palette = "Set2") +
  labs(fill = NULL) +
  annotate("text", x = venn_numbers$x, y = venn_numbers$y, label = venn_numbers$labels, size = 5) +
  annotate("text", x = 0, y = 2.65, label = "m4",fontface = "bold", size=4) +
  annotate("text", x = -2.3, y = -1.5,  label = "JC",fontface = "bold", size=4) +
  annotate("text", x = 2.3, y = -1.5, label = "m5",fontface = "bold", size=4) -> disagree_m1_venn
save_plot("temp.jpg", disagree_m1_venn)


pandit_au_plot <- plot_grid(pandit_au_barplot, disagree_m1_venn, labels="auto", scale=c(0.95, 0.8))
save_plot(paste0(figure_directory,"pandit_au.pdf"), pandit_au_plot, base_width=8, base_height=3.5)


######################### Model matrices and schematic ############################

###################################################################################
#### Plot schematic for m1-m5
# simulation_rf_fit %>% 
#   filter(name == scheme_name1, model == "poisson", rep == scheme_rep, tree == scheme_tree) %>% 
#   pull(BIC) -> poisson_bic_1
# simulation_rf_fit %>% 
#   filter(name == scheme_name1, model == "GTR20", rep == scheme_rep, tree == scheme_tree) %>% 
#   pull(BIC) -> gtr20_bic_1

# 
# all_sel %>% 
#   dplyr::select(name, tree, model, bic, repl) %>%
#   filter(name == scheme_name2, tree == scheme_tree, repl == scheme_rep) -> data2
# simulation_rf_fit %>% 
#   filter(name == scheme_name2, model == "poisson", rep == scheme_rep, tree == scheme_tree) %>% 
#   pull(BIC) -> poisson_bic_2
# simulation_rf_fit %>% 
#   filter(name == scheme_name2, model == "GTR20", rep == scheme_rep, tree == scheme_tree) %>% 
#   pull(BIC) -> gtr20_bic_2

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
  #bind_rows(
  #        tibble(name = scheme_name1, tree = scheme_tree, repl = scheme_rep, modelm = 4.3, modelm2 = "JC", model_name = "JC", model_levels_matrix = "JC", bic = poisson_bic_1),
  #        tibble(name = scheme_name1, tree = scheme_tree, repl = scheme_rep, modelm = 4.6, modelm2 = "GTR", model_name = "GTR", model_levels_matrix = "GTR", bic = gtr20_bic_1)
  #    ) %>%
  distinct() %>%
  mutate(model_levels= factor(modelm2, levels=c("m1", "m2", "m3", "m4", "m5"))) -> chunks

# msel_simulation %>%
#   filter(name == scheme_name2, tree == scheme_tree, repl == scheme_rep) %>%
#   rowwise() %>%
#   mutate(modelm2 = paste0("m", modelm),
#     model_levels_matrix = paste0(modelm2, " (", model_name, ")")) %>%
#   bind_rows(
#           tibble(name = scheme_name2, tree = scheme_tree, repl = scheme_rep, modelm = 4.3, modelm2 = "JC", model_name = "JC", model_levels_matrix = "JC", bic = poisson_bic_2),
#           tibble(name = scheme_name2, tree = scheme_tree, repl = scheme_rep, modelm = 4.6, modelm2 = "GTR", model_name = "GTR", model_levels_matrix = "GTR", bic = gtr20_bic_2)
#       ) %>%
#   distinct() %>%
#   mutate(model_levels= factor(modelm2, levels=model_labels)) -> chunks2

ggplot(data, aes(x = "", y = bic)) + 
    geom_violin(fill="dodgerblue3", color = "dodgerblue4", alpha=0.3) + 
    geom_point(alpha = 0.1, size=1.5) +
    geom_point(data = chunks, shape=21, aes(x = 1, y = bic, fill = model_levels), size=4) + 
    geom_label(data = chunks,  x = 1.02, hjust="outward", color = "black", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=4, fontface = "bold", alpha=0.8) + 
    #geom_label(data = subset(chunks, model_levels != "M1"),  x = 1.02, hjust="outward", color = "black", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.75) + 
    #geom_label(data = subset(chunks, model_levels == "M1"), x = 1.02, hjust="outward", color = "white", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.75) + 
    #ggtitle("HA simulation replicate") + 
    scale_fill_manual(values = c(  lighten(lighten(model_colors[1])), model_colors[2:7])) +
    theme(legend.position = "none", axis.ticks.x = element_blank()) +
    xlab("") + ylab("BIC values across all tested models")  -> quant_scheme_plot
ggsave(paste0(figure_directory, "bic_dist_ha_dosreis_ONLYm1-5.pdf"), quant_scheme_plot, width = 3.5, height=4)

# ggplot(data2, aes(x = "", y = bic)) + 
#     geom_violin(fill="dodgerblue3", color = "dodgerblue4", alpha=0.3) + 
#     geom_point(alpha = 0.1, size=1.5) +
#     geom_point(data = chunks2, shape=21, aes(x = 1, y = bic, fill = model_levels), size=3) + 
#     geom_label(data = chunks2,  x = 1.02, hjust="outward", color = "black", aes(y = bic + 25, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.8) + 
#     #geom_label(data = subset(chunks, model_levels != "M1"),  x = 1.02, hjust="outward", color = "black", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.75) + 
#     #geom_label(data = subset(chunks, model_levels == "M1"), x = 1.02, hjust="outward", color = "white", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.75) + 
#     ggtitle("1IBS simulation replicate") + 
#     scale_fill_manual(values = c(  lighten(lighten(model_colors[1])), model_colors[2:7])) +
#     theme(legend.position = "none", axis.ticks.x = element_blank()) +
#     xlab("") + ylab("BIC values across all tested models")  -> quant_scheme_plot_2

#bic_dist_grid <- plot_grid(quant_scheme_plot_1, quant_scheme_plot_2, labels="auto", nrow=1, scale=0.97)

#ggsave(paste0(figure_directory, "bic_dist_grid.pdf"), bic_dist_grid, width = 8, height=3)
    

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
    facet_wrap(~name_levels, nrow=1) + scale_fill_brewer(palette = "Dark2", name = paste("m1 Model Matrix")) +
    xlab("Simulation tree") + 
    ylab("Count") + 
    theme(legend.position = "bottom", axis.text.x = element_text(size=7, angle=20, margin = margin(t = 8, b=5))) -> selected_m1_plot        
save_plot(paste0(figure_directory, "selected_m1_simulation_barplot.pdf"), selected_m1_plot, base_width=10, base_height = 3.5) 


for (m in 2:5) {
    mname <- paste0("m", m)
    msel_simulation %>%
        filter(modelm == m) %>%
        rowwise() %>%
        mutate(model_matrix = str_split(model_name, "\\+")[[1]][1]) %>%
        group_by(name, tree, model_matrix) %>%
        tally() %>% 
        mutate(tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_abbr),
               name_levels  = factor(name, levels=name_levels),
               model_matrix_levels = fct_reorder(model_matrix, n)) %>% 
        ggplot(aes(x = tree_levels, y = n, fill = model_matrix_levels)) + 
        geom_bar(stat="identity") +
        facet_wrap(~name_levels, nrow=2) + 
        xlab("Simulation tree") + 
        ylab("Count") + 
        theme(axis.text.x = element_text(size=8, angle=30, margin = margin(t = 10, b = 5)), legend.position = "bottom")-> selected_m_plot        
        selected_m_plot <- selected_m_plot + guides(fill = guide_legend(nrow=2, title = paste(mname, "Model Matrix")))
            
    save_plot(paste0(figure_directory, "selected_", mname, "_simulation_barplot.pdf"), selected_m_plot, base_width=8, base_height = 4) 
}

      
    
x <- 1
plot_list <- c()
for (m in 1:5) {
    mname <- paste0("M", m)
    msel_pandit %>%
        filter(modelm == m) %>%
        rowwise() %>%
        mutate(model_matrix = str_split(model_name, "\\+")[[1]][1]) %>%
        group_by(model_matrix) %>%
        tally() %>% 
        mutate(model_matrix_levels = fct_reorder(model_matrix, n)) %>%
        ggplot(aes(x = model_matrix_levels, y = n, fill = model_matrix_levels)) + 
        geom_bar(stat="identity") +
        geom_text(aes(x = model_matrix_levels, y = n+2, label = n), size=3)+
        xlab("") + 
        ggtitle(paste(mname, "Model Matrix")) + 
        ylab("Count") + 
        #labs(fill = paste(mname, "Model Matrix")) + 
        theme(axis.text.x = element_text(size=8, angle=30, vjust=0.75),legend.position = "none") -> selected_m_plot
    if (m == 1) {
        selected_m_plot <- selected_m_plot + ggtitle("") + xlab(paste(mname, "Model Matrix")) 
        save_plot(paste0(figure_directory, "selected_", mname, "_pandit_barplot.pdf"), selected_m_plot, base_width=5, base_height = 3.25) 
    }
    if (m > 1) {
        plot_list[[x]] <- selected_m_plot
        x <- x+1
    } 
}
selected_grid <- plot_grid(plotlist = plot_list, nrow=2, labels = "auto")
save_plot(paste0(figure_directory, "selected_models_m2-5_pandit.pdf"), selected_grid, base_width=11, base_height=6) 



### entropy comparisons to distinguish simulation conditions
ggplot(entropy, aes(x = name_levels, y = entropy, color = name_levels)) + 
    geom_sina(alpha=0.7) + 
    scale_color_brewer(palette = "Set2") + 
    stat_summary(fun.data = "mean_se",size = 0.5, color = "black") + # se is too small to see so this ends up being mean, chachiiiiing 
    xlab("DMS simulation") + 
    ylab("Sitewise entropy") + 
    background_grid() + theme(legend.position = "none") -> entropy_sina
ggsave(paste0(figure_directory, "entropy_sina.pdf"), entropy_sina, width = 6, height=2.5)






