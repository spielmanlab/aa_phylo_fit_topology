source("load.R") # ugh i dont want to keep copy/pasting this stuff it lives there now
#stop()

figure_directory <- "evol2019_figures/"

######## simulation RF boxplots
simulation_rf_fit %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
         tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_ntaxa),
         name_levels  = factor(name, levels=name_levels, labels = name_labels_nsites)) -> sim_rf
         
sim_rf %>% 
    filter(tree %in% c("dosreis", "andersen"), name == "HA") %>%
        ggplot(aes(x = model_levels, y = rf_true_norm,fill = model_levels)) + 
            geom_boxplot(alpha=0, size=0.5) + 
            geom_point(position = position_jitterdodge(), shape=21)+ #geom_boxplot(outlier.size = 0.3, size=0.25,  width=0.9) + 
            scale_fill_manual(values = model_colors) +
            facet_wrap(~tree_levels, scales="free_y") + 
            panel_border() +
            background_grid() + 
            theme(legend.position = "none", panel.spacing = unit(20, "pt")) + 
            xlab("Protein Models") + ylab("Normalized RF distance") -> simulation_rf_boxplot  
save_plot(paste0(figure_directory,"simulation_rf_boxplot.pdf"), simulation_rf_boxplot, base_width=7, base_height=2.5)



sim_ufb %>% 
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
        tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_ntaxa),
        name_levels  = factor(name, levels=name_levels), 
        rep = factor(rep)) -> ufb_fact
ufb_fact %>%
    group_by(model_levels, tree_levels, name_levels, rep) %>%
    tally() %>%
    rename(total_considered = n) -> ufb_total_nodes
     
############## False positive nodes ##############
ufb_fact %>% 
    mutate(supported = boot >= 95) %>% 
    filter(in_true == FALSE, supported == TRUE) %>% 
    group_by(model_levels, tree_levels, name_levels, rep) %>%
    tally() %>% 
    complete(model_levels, tree_levels, name_levels, rep, fill = list(n = 0)) %>%
    left_join(ufb_total_nodes) %>%
    mutate(percent = n / total_considered) %>%
    distinct() -> ufb_fp

ufb_fact %>% 
    mutate(supported = boot >= 95) %>% 
    filter(in_true == supported) %>% 
    group_by(model_levels, tree_levels, name_levels, rep) %>%
    tally() %>% 
    complete(model_levels, tree_levels, name_levels, rep, fill = list(n = 0)) %>%
    left_join(ufb_total_nodes) %>%
    mutate(percent = n / total_considered) %>%
    distinct() -> ufb_accuracy



ufb_fp %>%
    filter(name_levels == "HA", tree_levels %in% c("Mammals (274)", "Lassa Virus (179)")) %>%
    ggplot(aes(x = model_levels, y = percent, fill = model_levels)) + 
        geom_boxplot(alpha=0, size=0.5) + 
        geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5,  size=2, color="grey10")  +     
        scale_fill_manual(values=model_colors, name = "Protein Model") +         
        facet_wrap(~tree_levels, nrow=1) +
        scale_y_continuous(breaks=c(0, 0.02, 0.04, 0.06)) + 
        background_grid() + 
        panel_border() +
        theme(legend.position = "none",
              panel.spacing = unit(20, "pt")) +
        xlab("Protein Models") + ylab("Proportion false positive nodes") +
        geom_hline(yintercept = 0.05, color = "dodgerblue3")  -> sim_ufb_fp_ha


ufb_accuracy %>%
    filter(name_levels == "HA", tree_levels %in% c("Mammals (274)", "Lassa Virus (179)")) %>%
    ggplot(aes(x = model_levels, y = percent, fill = model_levels)) + 
        geom_boxplot(alpha=0, size=0.5) + 
        geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5,  size=2, color="grey10")  +     
        scale_fill_manual(values=model_colors, name = "Protein Model") +         
        facet_wrap(~tree_levels, nrow=1, scales="free_y") +
        background_grid() + 
        panel_border() +
        theme(legend.position = "none",
              panel.spacing = unit(20, "pt")) +
        xlab("Protein Models") + ylab("Proportion accurate nodes") -> sim_ufb_acc_ha

ufb_fact %>%
    filter(rep == 1, name_levels == "HA", tree_levels %in% c("Mammals (274)", "Lassa Virus (179)")) %>%
    mutate(in_true_fct = factor(in_true, levels=c("TRUE", "FALSE"), labels=c("Yes", "No"))) %>%
    ggplot(aes(x = model_levels, y = boot, color = in_true_fct)) + 
        geom_boxplot(outlier.shape="", alpha=0) + 
        geom_point(position=position_jitterdodge(jitter.height=0, jitter.width=0.2),
                  alpha=0.5,  size=1)  +     
        #scale_fill_hue(l=50, name = "Node present in true tree")+
        scale_color_hue(l=50, name = "Node present in true tree")+
        facet_wrap(~tree_levels, nrow=1, scales="free_y") +
        background_grid() + 
        geom_hline(yintercept=95, color="black")+
        panel_border() +
        theme(legend.position = "bottom",
              panel.spacing = unit(20, "pt")) +
        xlab("Protein Models") + ylab("UFBoot2 Values") -> rep1_ufbbootdist
        


save_plot(paste0(figure_directory,"sim_ufb_fp_ha.pdf"), sim_ufb_fp_ha, base_width=9, base_height=3)
save_plot(paste0(figure_directory,"sim_ufb_acc_ha.pdf"), sim_ufb_acc_ha, base_width=9, base_height=3)
save_plot(paste0(figure_directory,"sim_rep1_ufbbootdistufb_fp_ha.pdf"), rep1_ufbbootdist, base_width=9, base_height=3)

        
        
        

#### Plot schematic for m1-m5
scheme_name <- "HA"
scheme_tree <- "ruhfel"
scheme_rep <- 1
all_sel %>% 
  dplyr::select(name, tree, model, bic, repl) %>%
  filter(name == scheme_name, tree == scheme_tree, repl == scheme_rep) -> scheme_data

simulation_rf_fit %>% 
  filter(name == scheme_name, model == "poisson", rep == scheme_rep, tree == scheme_tree) %>% 
  pull(BIC) -> poisson_bic

simulation_rf_fit %>% 
  filter(name == scheme_name, model == "GTR20", rep == scheme_rep, tree == scheme_tree) %>% 
  pull(BIC) -> gtr20_bic


msel_simulation %>%
  filter(name == scheme_name, tree == scheme_tree, repl == scheme_rep) %>%
  rowwise() %>%
  mutate(modelm2 = paste0("M", modelm),
    model_levels_matrix = paste0(modelm2, " (", model_name, ")")) %>%
  bind_rows(
          tibble(name = scheme_name, tree = scheme_tree, repl = scheme_rep, modelm = 4.3, modelm2 = "JC", model_name = "JC", model_levels_matrix = "JC", bic = poisson_bic),
          tibble(name = scheme_name, tree = scheme_tree, repl = scheme_rep, modelm = 4.6, modelm2 = "GTR", model_name = "GTR", model_levels_matrix = "GTR", bic = gtr20_bic)
      ) %>%
  mutate(model_levels= factor(modelm2, levels=c("M1", "M2", "M3", "M4", "M5", "JC", "GTR")))   -> chunks

chunks1 <- chunks %>% filter(!(model_levels %in% c("JC", "GTR")))
ggplot(scheme_data, aes(x = "", y = bic)) + 
    geom_violin(fill="dodgerblue3", color = "dodgerblue4", alpha=0.3) + 
  geom_point(alpha = 0.3, size=2.5) +
  geom_point(data = chunks1, shape=21, aes(x = 1, y = bic, fill = model_levels), size=4.5) + 
    geom_label(data = chunks1,  x = 1.02, hjust="outward", color = "black", 
                aes(y = bic+300, fill = model_levels, label = model_levels_matrix), 
                size=5, fontface = "bold", alpha=0.8) + 
    scale_y_continuous(limits=c(168000, 210000))+
    #geom_label(data = subset(chunks, model_levels != "M1"),  x = 1.02, hjust="outward", color = "black", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.75) + 
    #geom_label(data = subset(chunks, model_levels == "M1"), x = 1.02, hjust="outward", color = "white", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.75) + 
  scale_fill_manual(values = c(  lighten(lighten(model_colors[1])), model_colors[2:7])) +
  theme(axis.title = element_text(size=15), axis.text = element_text(size=14), legend.position = "none", axis.ticks.x = element_blank()) +
    xlab("") + ylab("BIC values across all tested models")  -> quant_scheme_plot

ggsave(paste0(figure_directory, "bic_dist_scheme_qviolin_m15.pdf"), quant_scheme_plot, width = 6, height=5)
            
        
        
stop()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
  # ufb_fact %>%
#     filter(rep == 1, name_levels == "HA", tree_levels %in% c("Mammals (274)", "Lassa Virus (179)")) %>%
#     mutate(in_true_fct = factor(in_true, levels=c("TRUE", "FALSE"), labels=c("Yes", "No"))) %>%
#     ggplot(aes(x = model_levels, y = boot, fill = in_true_fct, color = in_true_fct)) + 
#         geom_boxplot(outlier.shape="", alpha=0) + 
#         geom_point(position=position_jitterdodge(jitter.height=0, jitter.width=0.2), shape = 21, alpha=0.5,  size=1, color="grey10")  +     
#         scale_fill_hue(l=50, name = "Node present in true tree")+
#         scale_color_hue(l=50, name = "Node present in true tree")+
#         facet_wrap(~tree_levels, nrow=1, scales="free_y") +
#         background_grid() + 
#         geom_hline(yintercept=95)+
#         panel_border() +
#         theme(legend.position = "bottom",
#               panel.spacing = unit(20, "pt")) +
#         xlab("Protein Models") + ylab("UFBoot2 Values") -> rep1_ufbbootdist

    

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
  xlab("M1 Model Matrix") + 
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
    mutate(hasm1 = (model1 == "M1" | model2 == "M1")) %>%
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
        
        
    
    
    # theme(axis.line.y = element_blank(), 
    #           axis.ticks.y = element_blank(),
    #           legend.position = c(0.5, 0.01),  
    #           legend.key.height=unit(.7,"line"),
    #           legend.title = element_text(size=10), 
    #           legend.text = element_text(size=9), 
    #           axis.text = element_text(size=11),
    #           axis.title = element_text(size=12),
    #           axis.text.y = element_text(face = font_face) ) +

                                     


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
        scale_fill_hue(l=50, name = "In M1 Confidence Set") +
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
                      labels = c('M4', 'M5', 'JC'))
ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = 0.5, size = 1, colour = 'grey50') +
  theme_void() + 
  coord_fixed(clip = "off") + 
theme(legend.position = 'none') +
  #scale_fill_manual(values = c(poisson_col, m1to5_cols[4], m1to5_cols[5])) +
  #scale_colour_manual(values = c('cornflowerblue', 'firebrick', 'gold'), guide = FALSE) +
  scale_fill_hue(l=50) + labs(fill = NULL) +
  annotate("text", x = venn_numbers$x, y = venn_numbers$y, label = venn_numbers$labels, size = 5) +
  annotate("text", x = 0, y = 2.65, label = "M4",fontface = "bold", size=4) +
  annotate("text", x = -2.3, y = -1.5,  label = "JC",fontface = "bold", size=4) +
  annotate("text", x = 2.3, y = -1.5, label = "M5",fontface = "bold", size=4) -> disagree_m1_venn



pandit_au_plot <- plot_grid(pandit_au_barplot, disagree_m1_venn, labels="auto", scale=c(0.95, 0.8))
save_plot(paste0(figure_directory,"pandit_au.pdf"), pandit_au_plot, base_width=8, base_height=3.5)


#pandit_au_barplot_nolegend <- pandit_au_barplot + theme(legend.position = "none")
#pandit_au_barplot_legend<- get_legend(pandit_au_barplot)


#pandit_rf_nolegend <- rf_pandit_plot + theme(legend.position = "none")
#pandit_rf_legend<- get_legend(rf_pandit_plot)

#pandit_rf_au_plot <- plot_grid(pandit_rf_nolegend, pandit_au_barplot_nolegend, disagree_m1_venn, scale=c(1, 0.92, 0.7), nrow=1, labels="auto")

#part1 <- plot_grid(pandit_rf_nolegend, pandit_au_barplot_nolegend, nrow=1, labels=c("a", "b"))
#part2 <- plot_grid(pandit_rf_legend, pandit_au_barplot_legend, nrow=1)
#part12 <- plot_grid(part1, part2, nrow=2, rel_heights=c(1,0.1))
#pandit_rf_au_plot_full <- plot_grid(part12, disagree_m1_venn, nrow=1, labels=c("","c"),vjust=2.45, scale=c(0.95, 0.95, 0.5), rel_widths=c(0.7, 0.3))
#save_plot(paste0(figure_directory,"rf_pandit_plot.pdf"), pandit_rf_au_plot_full, base_width=14, base_height=5)


