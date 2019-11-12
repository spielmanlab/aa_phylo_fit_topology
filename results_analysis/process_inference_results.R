source("load.R") # ugh i dont want to keep copy/pasting this stuff it lives there now
#stop()


############################### Simulation figures ################################

#### simulation fit results
simulation_rf_fit %>%
  dplyr::select(-AIC, -AICc, -k,-logl, -rf, -rf_true_norm, -treelength, -max_rf, -ntaxa) %>%
  distinct() %>%
  group_by(name, tree, rep) %>%
  mutate(ic.rank = as.integer(rank(BIC))) %>%
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
           tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels),
           name_levels  = factor(name, levels=name_levels)) %>%
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
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
         tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_ntaxa),
         name_levels  = factor(name, levels=name_levels, labels = name_labels_nsites)) %>%
        ggplot(aes(x = model_levels, y = rf_true_norm)) + 
            #geom_sina(aes(color = name_levels), position = position_dodge(width=0.5)) + 
            geom_boxplot(aes(fill = name_levels), outlier.size = 0.3, size=0.25,  width=0.9) + 
            #geom_violin(size=0.15,  aes(fill = name_levels), scale = "width") + 
            #geom_point(position=position_jitterdodge(jitter.width = 0.2), size=0.75, alpha=0.7, aes(color = name_levels))  +       
            facet_wrap(~tree_levels, scales = "free_y", nrow=2) +
            #scale_color_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels) +
            #scale_fill_brewer(palette = "Set1", name = "DMS Simulations", labels = name_labels_nsites) +
            scale_fill_hue(l=50, name = "DMS Parameterization", labels = name_labels_nsites) +
            panel_border() +
            background_grid() + 
            #scale_y_continuous(limits = y_limits) + 
            xlab("Protein model") + ylab("Normalized Robinson-Foulds Distance") + 
            #geom_point(data = sim_rf_plotdata_mean, aes(x = model_levels, y = meanrf)) + 
            #geom_line(data = sim_rf_plotdata_mean, aes(x = model_levels, y = meanrf, group = name_levels)) + 
            theme(legend.position = "bottom") -> simulation_rf_boxplot  
save_plot(paste0(figure_directory,"simulation_rf_boxplot.pdf"), simulation_rf_boxplot, base_width=10, base_height=4)



#### NO YEAST GOT IT RIGHT ####
# simulation_rf_fit %>% group_by(name, tree, model) %>% tally(rf ==0) %>% filter(n>1) %>% xtable()
# % latex table generated in R 3.5.1 by xtable 1.8-4 package
# % Thu Jun  6 12:09:44 2019
# \begin{table}[ht]
# \centering
# \begin{tabular}{rlllr}
#   \hline
#  & name & tree & model & n \\ 
#   \hline
# 1 & HA & opisthokonta & GTR20 &   4 \\ 
#   2 & HA & salichos & GTR20 &   4 \\ 
#   3 & HA & salichos & m1 &   2 \\ 
#   4 & HA & salichos & m2 &   2 \\ 
#   5 & HIV & opisthokonta & GTR20 &  16 \\ 
#   6 & HIV & opisthokonta & m1 &   9 \\ 
#   7 & HIV & opisthokonta & m2 &   9 \\ 
#   8 & HIV & opisthokonta & m3 &   4 \\ 
#   9 & HIV & opisthokonta & m4 &   4 \\ 
#   10 & HIV & opisthokonta & poisson &   6 \\ 
#   11 & HIV & salichos & GTR20 &   6 \\ 
#   12 & HIV & salichos & m1 &   2 \\ 
#   13 & HIV & salichos & m2 &   2 \\ 
#   14 & HIV & salichos & m3 &   2 \\ 
#   15 & HIV & salichos & m4 &   2 \\ 
#   16 & HIV & salichos & poisson &   3 \\ 
#   17 & HIV & spiralia & GTR20 &   5 \\ 
#   18 & HIV & spiralia & m1 &   3 \\ 
#   19 & HIV & spiralia & m2 &   2 \\ 
#   20 & HIV & spiralia & m3 &   3 \\ 
#   21 & HIV & spiralia & m4 &   4 \\ 
#   22 & HIV & spiralia & poisson &   4 \\ 
#   23 & NP & opisthokonta & GTR20 &   3 \\ 
#   24 & NP & salichos & GTR20 &   4 \\ 
#   25 & NP & salichos & m1 &   4 \\ 
#   26 & NP & salichos & m2 &   3 \\ 
#   27 & NP & salichos & m4 &   3 \\     
# \hline
# \end{tabular}
# \end{table}



simulation_topology %>%
    filter(whichtest == "au") %>%
    gather(model, pvalue, m1:true) %>%
    mutate(sig = pvalue < 0.01) %>%
    filter(sig == TRUE) %>%
    group_by(name, tree, model) %>%
    tally() %>%
    group_by(model) %>% tally()
#1 GTR20     1
#2 m3        1
#3 m4        2
#4 m5       10
#5 true     29

simulation_topology %>%
    filter(whichtest == "sh") %>%
    gather(model, pvalue, m1:true) %>%
    mutate(sig = pvalue < 0.01) %>%
    filter(sig == TRUE) %>%
    group_by(name, tree, model) %>%
    tally() %>%
#1 m5        2
#2 true      9



     
############## False positive nodes ##############
ufb_fact %>% 
    count(model_levels, tree_levels, name_levels, rep, classif) %>%
    #complete(model_levels, tree_levels, name_levels, rep, classif, fill = list(n = 0)) %>%
    pivot_wider(names_from = classif, values_from = n) %>%
    ungroup() %>%
    replace_na(list(FP = 0, FN = 0, TP = 0, TN = 0)) %>%
    mutate(FPR = ifelse(is.nan(FP / (TN+FP)), 0, FP / (TN+FP)), 
           accuracy = (TP + TN)/(TP+TN+FP+FN)) -> ufb_fact_classif

some_trees <- c("Green Plant (360)", "Aves (200)","Opisthokonta (70)")
some_names <- c("1IBS", "HA")
ufb_fact_classif %>%
    filter(name_levels %in% some_names, tree_levels %in% some_trees) %>%
    ggplot(aes(x = model_levels, y = FPR, fill = model_levels)) + 
        geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5, color="grey10")  +     
        scale_fill_manual(values=model_colors, name = "Protein Model") +         
        facet_grid(name_levels~tree_levels, scales="free_y") +
        #scale_y_continuous(breaks=c(0, 0.02, 0.04, 0.06)) + 
        background_grid() + 
        panel_border() +
        theme(legend.position = "none",
              strip.text = element_text(size=8), 
              strip.background = element_rect(fill = "grey80"), 
              axis.text.x = element_text(size=8),
              panel.spacing = unit(0.2, "cm")) +
        xlab("Protein Models") + ylab("False positive rate") +
        geom_hline(yintercept = 0.05, color = "dodgerblue3")  -> fpr_subset

ufb_fact_classif %>%
    filter(name_levels %in% some_names, tree_levels %in% some_trees) %>%
    ggplot(aes(x = model_levels, y = accuracy, fill = model_levels)) + 
        geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5, color="grey10")  +     
        scale_fill_manual(values=model_colors, name = "Protein Model") +         
        facet_grid(name_levels~tree_levels, scales="free_y") +
        #scale_y_continuous(breaks=c(0, 0.02, 0.04, 0.06)) + 
        background_grid() + 
        panel_border() +
        theme(legend.position = "none",
              strip.text = element_text(size=8), 
              strip.background = element_rect(fill = "grey80"), 
              axis.text.x = element_text(size=8),
              panel.spacing = unit(0.2, "cm")) +
        xlab("Protein Models") + ylab("Accuracy")   -> acc_subset


ufb_subset_grid <- plot_grid(fpr_subset, acc_subset, nrow=2, labels= "auto", scale=0.98)
save_plot(paste0(figure_directory,"ufb_subset_grid.pdf"), ufb_subset_grid, base_width = 6, base_height=7)
             
             
             
             
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
              strip.background = element_rect(fill = "grey80"), 
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
              strip.background = element_rect(fill = "grey80"), 
              axis.text.x = element_text(size=8),
              panel.spacing = unit(0.2, "cm")) +
        xlab("Protein Models") + ylab("Accuracy")   -> acc_all
save_plot(paste0(figure_directory,"ufb_fpr_all.pdf"), fpr_all, base_width = 12, base_height=10)
save_plot(paste0(figure_directory,"ufb_acc_all.pdf"), acc_all, base_width = 12, base_height=10)


######### TODO: THE SI VERSION WITH ALL DATA ########


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
    mutate(notsig = pvalue >= 0.01) %>%
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
        

pandit_topology %>% filter(whichtest == "au", m4 >=0.01, m5>=0.01, poisson <0.01) %>% nrow()  -> poisson_only   
pandit_topology %>% filter(whichtest == "au", m4 >=0.01, m5<0.01, poisson >=0.01) %>% nrow()  -> m5_only
pandit_topology %>% filter(whichtest == "au", m4 <0.01,  m5>=0.01, poisson >=0.01) %>% nrow() -> m4_only
pandit_topology %>% filter(whichtest == "au", m4 >=0.01, m5<0.01, poisson <0.01) %>% nrow()   -> m5_poisson       
pandit_topology %>% filter(whichtest == "au", m4 <0.01,  m5<0.01, poisson >=0.01) %>% nrow()  -> m4_m5
pandit_topology %>% filter(whichtest == "au", m4 <0.01,  m5>=0.01, poisson <0.01) %>% nrow()  -> m4_poisson
pandit_topology %>% filter(whichtest == "au", m4 <0.01,  m5<0.01, poisson <0.01) %>% nrow()   -> all_three        
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
scheme_name1 <- "HA"
scheme_name2 <- "1IBS"
scheme_tree <- "dosreis"
scheme_rep <- 7
all_sel %>% 
  dplyr::select(name, tree, model, bic, repl) %>%
  filter(name == scheme_name1, tree == scheme_tree, repl == scheme_rep) -> data1
simulation_rf_fit %>% 
  filter(name == scheme_name1, model == "poisson", rep == scheme_rep, tree == scheme_tree) %>% 
  pull(BIC) -> poisson_bic_1
simulation_rf_fit %>% 
  filter(name == scheme_name1, model == "GTR20", rep == scheme_rep, tree == scheme_tree) %>% 
  pull(BIC) -> gtr20_bic_1


all_sel %>% 
  dplyr::select(name, tree, model, bic, repl) %>%
  filter(name == scheme_name2, tree == scheme_tree, repl == scheme_rep) -> data2
simulation_rf_fit %>% 
  filter(name == scheme_name2, model == "poisson", rep == scheme_rep, tree == scheme_tree) %>% 
  pull(BIC) -> poisson_bic_2
simulation_rf_fit %>% 
  filter(name == scheme_name2, model == "GTR20", rep == scheme_rep, tree == scheme_tree) %>% 
  pull(BIC) -> gtr20_bic_2


msel_simulation %>%
  filter(name == scheme_name1, tree == scheme_tree, repl == scheme_rep) %>%
  rowwise() %>%
  mutate(modelm2 = paste0("m", modelm),
    model_levels_matrix = paste0(modelm2, " (", model_name, ")")) %>%
  bind_rows(
          tibble(name = scheme_name1, tree = scheme_tree, repl = scheme_rep, modelm = 4.3, modelm2 = "JC", model_name = "JC", model_levels_matrix = "JC", bic = poisson_bic_1),
          tibble(name = scheme_name1, tree = scheme_tree, repl = scheme_rep, modelm = 4.6, modelm2 = "GTR", model_name = "GTR", model_levels_matrix = "GTR", bic = gtr20_bic_1)
      ) %>%
  distinct() %>%
  mutate(model_levels= factor(modelm2, levels=model_labels)) -> chunks1

msel_simulation %>%
  filter(name == scheme_name2, tree == scheme_tree, repl == scheme_rep) %>%
  rowwise() %>%
  mutate(modelm2 = paste0("m", modelm),
    model_levels_matrix = paste0(modelm2, " (", model_name, ")")) %>%
  bind_rows(
          tibble(name = scheme_name2, tree = scheme_tree, repl = scheme_rep, modelm = 4.3, modelm2 = "JC", model_name = "JC", model_levels_matrix = "JC", bic = poisson_bic_2),
          tibble(name = scheme_name2, tree = scheme_tree, repl = scheme_rep, modelm = 4.6, modelm2 = "GTR", model_name = "GTR", model_levels_matrix = "GTR", bic = gtr20_bic_2)
      ) %>%
  distinct() %>%
  mutate(model_levels= factor(modelm2, levels=model_labels)) -> chunks2

ggplot(data1, aes(x = "", y = bic)) + 
    geom_violin(fill="dodgerblue3", color = "dodgerblue4", alpha=0.3) + 
    geom_point(alpha = 0.1, size=1.5) +
    geom_point(data = chunks1, shape=21, aes(x = 1, y = bic, fill = model_levels), size=3) + 
    geom_label(data = chunks1,  x = 1.02, hjust="outward", color = "black", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.8) + 
    #geom_label(data = subset(chunks, model_levels != "M1"),  x = 1.02, hjust="outward", color = "black", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.75) + 
    #geom_label(data = subset(chunks, model_levels == "M1"), x = 1.02, hjust="outward", color = "white", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.75) + 
    ggtitle("HA simulation replicate") + 
    scale_fill_manual(values = c(  lighten(lighten(model_colors[1])), model_colors[2:7])) +
    theme(legend.position = "none", axis.ticks.x = element_blank()) +
    xlab("") + ylab("BIC values across all tested models")  -> quant_scheme_plot_1

ggplot(data2, aes(x = "", y = bic)) + 
    geom_violin(fill="dodgerblue3", color = "dodgerblue4", alpha=0.3) + 
    geom_point(alpha = 0.1, size=1.5) +
    geom_point(data = chunks2, shape=21, aes(x = 1, y = bic, fill = model_levels), size=3) + 
    geom_label(data = chunks2,  x = 1.02, hjust="outward", color = "black", aes(y = bic + 25, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.8) + 
    #geom_label(data = subset(chunks, model_levels != "M1"),  x = 1.02, hjust="outward", color = "black", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.75) + 
    #geom_label(data = subset(chunks, model_levels == "M1"), x = 1.02, hjust="outward", color = "white", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.75) + 
    ggtitle("1IBS simulation replicate") + 
    scale_fill_manual(values = c(  lighten(lighten(model_colors[1])), model_colors[2:7])) +
    theme(legend.position = "none", axis.ticks.x = element_blank()) +
    xlab("") + ylab("BIC values across all tested models")  -> quant_scheme_plot_2

bic_dist_grid <- plot_grid(quant_scheme_plot_1, quant_scheme_plot_2, labels="auto", nrow=1, scale=0.97)

ggsave(paste0(figure_directory, "bic_dist_grid.pdf"), bic_dist_grid, width = 8, height=3)
    
    
    
    
    
    

for (m in 1:5) {
    mname <- paste0("M", m)
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
        theme(axis.text.x = element_text(size=8, angle=30, margin = margins(t = 10, b = 5)), legend.position = "bottom")-> selected_m_plot        
        if (m == 1)
        {
            selected_m_plot <- selected_m_plot + scale_fill_hue(l=50, name = paste(mname, "Model Matrix")) 
        }   else{
            selected_m_plot <- selected_m_plot + guides(fill = guide_legend(nrow=2, title = paste(mname, "Model Matrix")))
        }
            
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
    xlab("Simulation parameterization") + 
    ylab("Sitewise entropy") + 
    background_grid() + theme(legend.position = "none") -> entropy_sina
ggsave(paste0(figure_directory, "entropy_sina.pdf"), entropy_sina, width = 6, height=2.5)

############################### Linear models ################################

######## RF for simulations
simulation_rf_fit$model <- factor(simulation_rf_fit$model, levels=c("m1", "m2", "m3", "m4", "m5", "poisson", "GTR20"))
lmer(rf_true_norm ~ model + (1|name) + (1|tree), data = simulation_rf_fit) -> fit
glht(fit, linfct=mcp(model='Tukey')) %>% summary()
#                        Estimate Std. Error z value Pr(>|z|)    
# m2 - m1 == 0          0.0028237  0.0034829   0.811  0.98398    
# m3 - m1 == 0          0.0050753  0.0034829   1.457  0.77026    
# m4 - m1 == 0          0.0077192  0.0034829   2.216  0.28652    
# poisson - m1 == 0     0.0026428  0.0034829   0.759  0.98867    
# GTR20 - m1 == 0      -0.0024479  0.0034829  -0.703  0.99247    
# m3 - m2 == 0          0.0022516  0.0034829   0.646  0.99522    
# m4 - m2 == 0          0.0048955  0.0034829   1.406  0.79917    
# poisson - m2 == 0    -0.0001809  0.0034829  -0.052  1.00000    
# GTR20 - m2 == 0      -0.0052716  0.0034829  -1.514  0.73699    
# m4 - m3 == 0          0.0026439  0.0034829   0.759  0.98865    
# poisson - m3 == 0    -0.0024325  0.0034829  -0.698  0.99272    
# GTR20 - m3 == 0      -0.0075232  0.0034829  -2.160  0.31765    
# poisson - m4 == 0    -0.0050764  0.0034829  -1.458  0.76998    
# GTR20 - m4 == 0      -0.0101671  0.0034829  -2.919  0.05429 .  
# GTR20 - poisson == 0 -0.0050907  0.0034829  -1.462  0.76767 

# m5 - m1 == 0          0.0216532  0.0034829   6.217  < 0.001 ***
# m5 - m2 == 0          0.0188295  0.0034829   5.406  < 0.001 ***
# m5 - m3 == 0          0.0165779  0.0034829   4.760  < 0.001 ***
# m5 - m4 == 0          0.0139340  0.0034829   4.001  0.00127 ** 
# poisson - m5 == 0    -0.0190104  0.0034829  -5.458  < 0.001 ***
# GTR20 - m5 == 0      -0.0241011  0.0034829  -6.920  < 0.001 ***




############################ linear model fp, acc ##################
#fit <- lmer(percent ~ model_levels + (1|tree_levels) + (1|name_levels), data = ufb_fp)
#glht(fit, linfct=mcp(model_levels='Tukey')) %>% summary()
# 
# 	 Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: lmer(formula = percent ~ model_levels + (1 | tree_levels) + (1 | 
#     name_levels), data = ufb_fp)
# 
# Linear Hypotheses:
#                 Estimate Std. Error z value Pr(>|z|)    
# m2 - m1 == 0   1.079e-03  6.815e-04   1.583  0.69375    
# m3 - m1 == 0   2.392e-03  6.815e-04   3.510  0.00815 ** 
# m4 - m1 == 0   3.851e-03  6.815e-04   5.650  < 0.001 ***
# m5 - m1 == 0   1.055e-02  6.815e-04  15.475  < 0.001 ***
# JC - m1 == 0   2.325e-03  6.815e-04   3.412  0.01149 *  NS
# GTR - m1 == 0  7.051e-03  6.815e-04  10.346  < 0.001 ***
# m3 - m2 == 0   1.314e-03  6.815e-04   1.928  0.46186    

# m4 - m2 == 0   2.772e-03  6.815e-04   4.068  < 0.001 ***
# m5 - m2 == 0   9.468e-03  6.815e-04  13.892  < 0.001 ***
# JC - m2 == 0   1.247e-03  6.815e-04   1.829  0.52812    
# GTR - m2 == 0  5.973e-03  6.815e-04   8.764  < 0.001 ***


# m4 - m3 == 0   1.458e-03  6.815e-04   2.140  0.32896    
# m5 - m3 == 0   8.154e-03  6.815e-04  11.964  < 0.001 ***
# JC - m3 == 0  -6.715e-05  6.815e-04  -0.099  1.00000    
# GTR - m3 == 0  4.659e-03  6.815e-04   6.836  < 0.001 ***

# m5 - m4 == 0   6.696e-03  6.815e-04   9.825  < 0.001 ***
# JC - m4 == 0  -1.526e-03  6.815e-04  -2.238  0.27471    
# GTR - m4 == 0  3.200e-03  6.815e-04   4.696  < 0.001 ***
# JC - m5 == 0  -8.221e-03  6.815e-04 -12.063  < 0.001 ***
# GTR - m5 == 0 -3.495e-03  6.815e-04  -5.129  < 0.001 ***
# GTR - JC == 0  4.726e-03  6.815e-04   6.934  < 0.001 ***


fit <- lmer(percent ~ model_levels + (1|tree_levels) + (1|name_levels), data = ufb_accuracy)
glht(fit, linfct=mcp(model_levels='Tukey')) %>% summary()
# Linear Hypotheses:
#                 Estimate Std. Error z value Pr(>|z|)    
# m2 - m1 == 0  -0.0009798  0.0023054  -0.425  0.99955    
# m3 - m1 == 0   0.0005854  0.0023054   0.254  0.99998    
# m4 - m1 == 0  -0.0015766  0.0023054  -0.684  0.99351    
# m5 - m1 == 0  -0.0048861  0.0023054  -2.119  0.34129    
# JC - m1 == 0  -0.0020356  0.0023054  -0.883  0.97520    
# GTR - m1 == 0  0.0084079  0.0023054   3.647  0.00505 ** 

# m3 - m2 == 0   0.0015652  0.0023054   0.679  0.99376    
# m4 - m2 == 0  -0.0005968  0.0023054  -0.259  0.99998    
# m5 - m2 == 0  -0.0039063  0.0023054  -1.694  0.61975    
# JC - m2 == 0  -0.0010558  0.0023054  -0.458  0.99931    
# GTR - m2 == 0  0.0093876  0.0023054   4.072  < 0.001 ***

# m4 - m3 == 0  -0.0021620  0.0023054  -0.938  0.96646    
# m5 - m3 == 0  -0.0054715  0.0023054  -2.373  0.20995    
# JC - m3 == 0  -0.0026210  0.0023054  -1.137  0.91677    
# GTR - m3 == 0  0.0078225  0.0023054   3.393  0.01222 * NS
 
# m5 - m4 == 0  -0.0033095  0.0023054  -1.436  0.78248    
# JC - m4 == 0  -0.0004590  0.0023054  -0.199  0.99999    
# GTR - m4 == 0  0.0099844  0.0023054   4.331  < 0.001 ***
# JC - m5 == 0   0.0028505  0.0023054   1.236  0.87999    
# GTR - m5 == 0  0.0132939  0.0023054   5.767  < 0.001 ***
# GTR - JC == 0  0.0104435  0.0023054   4.530  < 0.001 ***


# 

