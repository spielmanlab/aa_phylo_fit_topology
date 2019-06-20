source("load.R") # ugh i dont want to keep copy/pasting this stuff it lives there now
#stop()


############################### Simulation figures ################################

#### simulation fit results
simulation_rf_fit %>%
  dplyr::select(-AIC, -AICc, -k,-logl, -rf, -rf_true_norm, -treelength) %>%
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
    filter(name_levels == "HA") %>%
    ggplot(aes(x = model_levels, y = percent, fill = model_levels)) + 
        geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5,  size=2, color="grey10")  +     
        scale_fill_manual(values=model_colors, name = "Protein Model") +         
        facet_wrap(~tree_levels, nrow=2) +
        scale_y_continuous(breaks=c(0, 0.02, 0.04, 0.06)) + 
        background_grid() + 
        panel_border() +
        theme(legend.position = "none",
              strip.text = element_text(size=8), 
              axis.text.x = element_text(size=8)) +
        xlab("Protein Models") + ylab("Proportion false positive nodes") +
        geom_hline(yintercept = 0.05, color = "dodgerblue3")  -> sim_ufb_fp_ha


ufb_accuracy %>%
    filter(name_levels == "HA") %>%
    ggplot(aes(x = model_levels, y = percent, fill = model_levels)) + 
        geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5,  size=2, color="grey10")  +     
        scale_fill_manual(values=model_colors, name = "Protein Model") +         
        facet_wrap(~tree_levels, nrow=2, scales="free_y") +
        background_grid() + 
        panel_border() +
        theme(legend.position = "none",
              strip.text = element_text(size=8), 
              axis.text = element_text(size=8)) +
             xlab("Protein Models") + ylab("Proportion accurate nodes") -> sim_ufb_acc_ha


ufb_ha_grid <- plot_grid(sim_ufb_fp_ha, sim_ufb_acc_ha, nrow=2, labels="auto", scale=0.98)
save_plot(paste0(figure_directory,"ufb_simulation_fp_accuracy_HA.pdf"), ufb_ha_grid, base_width = 9, base_height=5)




ufb_fp %>%
    filter(name_levels %in% c("HIV", "NP")) %>%
    ggplot(aes(x = model_levels, y = percent, fill = model_levels)) + 
        geom_jitter(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5,  size=1.75, color="grey10")  +     
        scale_fill_manual(values=model_colors, name = "Protein Model") +         
        facet_grid(name_levels~tree_levels) +
        background_grid() + 
        panel_border() +
        theme(legend.position = "none",
              strip.text = element_text(size=8), 
              axis.text.x = element_text(size=6, angle=30)) +
        xlab("Protein Models") + ylab("Percent false positive nodes") +
        geom_hline(yintercept = 0.05, color = "dodgerblue3")  -> sim_ufb_fp_nphiv
save_plot(paste0(figure_directory,"ufb_simulation_fp_HIV-NP.pdf"), sim_ufb_fp_nphiv, base_width = 10, base_height=3)




ufb_accuracy %>%
    filter(name_levels %in% c("HIV", "NP")) %>%
    ggplot(aes(x = model_levels, y = percent, fill = model_levels)) + 
        geom_point(position=position_jitterdodge(jitter.height=0), shape = 21, alpha=0.5,  size=1.75, color="grey10")  +     
        scale_fill_manual(values=model_colors, name = "Protein Model") +         
        facet_grid(name_levels~tree_levels) +
        scale_y_continuous(breaks=seq(0.5, 1.0, 0.1)) + 
        background_grid() + 
        panel_border() +
        theme(legend.position = "none",
              strip.text = element_text(size=8), 
              axis.text.x = element_text(size=6, angle=30)) +
        xlab("Protein Models") + ylab("Percent accurate nodes") -> sim_ufb_acc_nphiv
save_plot(paste0(figure_directory,"ufb_simulation_accuracy_HIV-NP.pdf"), sim_ufb_acc_nphiv, base_width = 10, base_height=3)

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
# M2 - M1 == 0   0.0008606  0.0006003   1.434  0.78363    
# M3 - M1 == 0   0.0017267  0.0006003   2.876  0.06123 .  
# M4 - M1 == 0   0.0015283  0.0006003   2.546  0.14294    
# M5 - M1 == 0   0.0060008  0.0006003   9.996  < 0.001 ***
# JC - M1 == 0   0.0009377  0.0006003   1.562  0.70661    
# GTR - M1 == 0  0.0022484  0.0006003   3.745  0.00342 ** 
# M3 - M2 == 0   0.0008661  0.0006003   1.443  0.77846    
# M4 - M2 == 0   0.0006676  0.0006003   1.112  0.92462    
# M5 - M2 == 0   0.0051402  0.0006003   8.563  < 0.001 ***
# JC - M2 == 0   0.0000771  0.0006003   0.128  1.00000    
# GTR - M2 == 0  0.0013878  0.0006003   2.312  0.23789    
# M4 - M3 == 0  -0.0001985  0.0006003  -0.331  0.99990    
# M5 - M3 == 0   0.0042741  0.0006003   7.120  < 0.001 ***
# JC - M3 == 0  -0.0007891  0.0006003  -1.314  0.84554    
# GTR - M3 == 0  0.0005216  0.0006003   0.869  0.97710    
# M5 - M4 == 0   0.0044726  0.0006003   7.451  < 0.001 ***
# JC - M4 == 0  -0.0005906  0.0006003  -0.984  0.95763    
# GTR - M4 == 0  0.0007201  0.0006003   1.200  0.89459    
# JC - M5 == 0  -0.0050631  0.0006003  -8.434  < 0.001 ***
# GTR - M5 == 0 -0.0037525  0.0006003  -6.251  < 0.001 ***
# GTR - JC == 0  0.0013107  0.0006003   2.183  0.30437  

# M5 - M1 == 0   0.0060008  0.0006003   9.996  < 0.001 ***
# GTR - M1 == 0  0.0022484  0.0006003   3.745  0.00342 ** 
# M5 - M2 == 0   0.0051402  0.0006003   8.563  < 0.001 ***
# JC - M5 == 0  -0.0050631  0.0006003  -8.434  < 0.001 ***
# M5 - M4 == 0   0.0044726  0.0006003   7.451  < 0.001 ***
# M5 - M3 == 0   0.0042741  0.0006003   7.120  < 0.001 ***
# GTR - M5 == 0 -0.0037525  0.0006003  -6.251  < 0.001 ***


fit <- lmer(percent ~ model_levels + (1|tree_levels) + (1|name_levels), data = ufb_accuracy)
glht(fit, linfct=mcp(model_levels='Tukey')) %>% summary()
#                 Estimate Std. Error z value Pr(>|z|)    
# M2 - M1 == 0  -0.0022988  0.0024671  -0.932  0.96752    
# M3 - M1 == 0  -0.0038431  0.0024671  -1.558  0.70950    
# M4 - M1 == 0  -0.0056102  0.0024671  -2.274  0.25679    
# M5 - M1 == 0  -0.0156591  0.0024671  -6.347  < 0.001 ***
# JC - M1 == 0  -0.0062340  0.0024671  -2.527  0.14969    
# GTR - M1 == 0  0.0099124  0.0024671   4.018  0.00125 ** 
# M3 - M2 == 0  -0.0015443  0.0024671  -0.626  0.99600    
# M4 - M2 == 0  -0.0033115  0.0024671  -1.342  0.83195    
# M5 - M2 == 0  -0.0133603  0.0024671  -5.415  < 0.001 ***
# JC - M2 == 0  -0.0039352  0.0024671  -1.595  0.68552    
# GTR - M2 == 0  0.0122111  0.0024671   4.950  < 0.001 ***
# M4 - M3 == 0  -0.0017671  0.0024671  -0.716  0.99167    
# M5 - M3 == 0  -0.0118159  0.0024671  -4.789  < 0.001 ***
# JC - M3 == 0  -0.0023908  0.0024671  -0.969  0.96061    
# GTR - M3 == 0  0.0137555  0.0024671   5.575  < 0.001 ***
# M5 - M4 == 0  -0.0100488  0.0024671  -4.073  < 0.001 ***
# JC - M4 == 0  -0.0006237  0.0024671  -0.253  0.99998    
# GTR - M4 == 0  0.0155226  0.0024671   6.292  < 0.001 ***
# JC - M5 == 0   0.0094251  0.0024671   3.820  0.00251 ** 
# GTR - M5 == 0  0.0255714  0.0024671  10.365  < 0.001 ***
# GTR - JC == 0  0.0161463  0.0024671   6.545  < 0.001 ***
# ---

# M5 - M1 == 0  -0.0156591  0.0024671  -6.347  < 0.001 ***
# M5 - M2 == 0  -0.0133603  0.0024671  -5.415  < 0.001 ***
# M5 - M3 == 0  -0.0118159  0.0024671  -4.789  < 0.001 ***
# M5 - M4 == 0  -0.0100488  0.0024671  -4.073  < 0.001 ***
# JC - M5 == 0   0.0094251  0.0024671   3.820  0.00251 ** 
# GTR - M1 == 0  0.0099124  0.0024671   4.018  0.00125 ** 
# GTR - M2 == 0  0.0122111  0.0024671   4.950  < 0.001 ***
# GTR - M3 == 0  0.0137555  0.0024671   5.575  < 0.001 ***
# GTR - M4 == 0  0.0155226  0.0024671   6.292  < 0.001 ***
# GTR - M5 == 0  0.0255714  0.0024671  10.365  < 0.001 ***
# GTR - JC == 0  0.0161463  0.0024671   6.545  < 0.001 ***


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
        xlab("GTR20 model rank") + ylab("Number of taxa") + 
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




######################### Model matrices and schematic ############################

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
          tibble(name = scheme_name, tree = scheme_tree, repl = scheme_rep, modelm = 4.6, modelm2 = "GTR20", model_name = "GTR20", model_levels_matrix = "GTR20", bic = gtr20_bic)
      ) %>%
  mutate(model_levels= factor(modelm2, levels=c("M1", "M2", "M3", "M4", "M5", "JC", "GTR20")))   -> chunks


ggplot(scheme_data, aes(x = "", y = bic)) + 
    geom_violin(fill="dodgerblue3", color = "dodgerblue4", alpha=0.3) + 
  geom_point(alpha = 0.3, size=2.5) +
  geom_point(data = chunks, shape=21, aes(x = 1, y = bic, fill = model_levels), size=2.5) + 
    geom_label(data = chunks,  x = 1.02, hjust="outward", color = "black", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.8) + 
    #geom_label(data = subset(chunks, model_levels != "M1"),  x = 1.02, hjust="outward", color = "black", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.75) + 
    #geom_label(data = subset(chunks, model_levels == "M1"), x = 1.02, hjust="outward", color = "white", aes(y = bic+300, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.75) + 
  scale_fill_manual(values = c(  lighten(lighten(model_colors[1])), model_colors[2:7])) +
  theme(legend.position = "none", axis.ticks.x = element_blank()) +
    xlab("") + ylab("BIC values across all tested models")  -> quant_scheme_plot

ggsave(paste0(figure_directory, "bic_dist_scheme_qviolin.pdf"), quant_scheme_plot, width = 6, height=5)
    

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
        facet_wrap(~name_levels, nrow=1) + 
        xlab("Simulation tree") + 
        ylab("Count") + 
        theme(axis.text.x = element_text(size=8), legend.position = "bottom")-> selected_m_plot        
        if (m == 1)
        {
            selected_m_plot <- selected_m_plot + scale_fill_hue(l=50, name = paste(mname, "Model Matrix")) 
        }   else{
            selected_m_plot <- selected_m_plot + guides(fill = guide_legend(nrow=2, title = paste(mname, "Model Matrix")))
        }
            
    save_plot(paste0(figure_directory, "selected_", mname, "_simulation_barplot.pdf"), selected_m_plot, base_width=10, base_height = 3) 
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



############################### Linear models ################################

######## RF for simulations
simulation_rf_fit$model <- factor(simulation_rf_fit$model, levels=c("m1", "m2", "m3", "m4", "m5", "poisson", "GTR20"))
lmer(rf_true_norm ~ model + (1|name) + (1|tree), data = simulation_rf_fit) -> fit
glht(fit, linfct=mcp(model='Tukey')) %>% summary()
# Linear Hypotheses:
# m2 - m1 == 0          0.0024572  0.0020367   1.206  0.89194    
# m3 - m1 == 0          0.0044539  0.0020367   2.187  0.30229    
# m4 - m1 == 0          0.0064229  0.0020367   3.154  0.02686 *  
# m5 - m1 == 0          0.0226019  0.0020367  11.097  < 0.001 ***
# poisson - m1 == 0     0.0034234  0.0020367   1.681  0.62890    
# GTR20 - m1 == 0      -0.0072462  0.0020367  -3.558  0.00693 ** 
# m3 - m2 == 0          0.0019968  0.0020367   0.980  0.95833    
# m4 - m2 == 0          0.0039657  0.0020367   1.947  0.44881    
# m5 - m2 == 0          0.0201447  0.0020367   9.891  < 0.001 ***
# poisson - m2 == 0     0.0009663  0.0020367   0.474  0.99916    
# GTR20 - m2 == 0      -0.0097033  0.0020367  -4.764  < 0.001 ***
# m4 - m3 == 0          0.0019690  0.0020367   0.967  0.96108    
# m5 - m3 == 0          0.0181480  0.0020367   8.910  < 0.001 ***
# poisson - m3 == 0    -0.0010305  0.0020367  -0.506  0.99878    
# GTR20 - m3 == 0      -0.0117001  0.0020367  -5.745  < 0.001 ***
# m5 - m4 == 0          0.0161790  0.0020367   7.944  < 0.001 ***
# poisson - m4 == 0    -0.0029994  0.0020367  -1.473  0.76118    
# GTR20 - m4 == 0      -0.0136690  0.0020367  -6.711  < 0.001 ***
# poisson - m5 == 0    -0.0191785  0.0020367  -9.416  < 0.001 ***
# GTR20 - m5 == 0      -0.0298481  0.0020367 -14.655  < 0.001 ***
# GTR20 - poisson == 0 -0.0106696  0.0020367  -5.239  < 0.001 ***

## SIG LINES:
# GTR20 - m1 == 0      -0.0072462  0.0020367  -3.558  0.00693 ** 
# GTR20 - m2 == 0      -0.0097033  0.0020367  -4.764  < 0.001 ***
# GTR20 - m3 == 0      -0.0117001  0.0020367  -5.745  < 0.001 ***
# GTR20 - m5 == 0      -0.0298481  0.0020367 -14.655  < 0.001 ***
# GTR20 - m4 == 0      -0.0136690  0.0020367  -6.711  < 0.001 ***
# GTR20 - poisson == 0 -0.0106696  0.0020367  -5.239  < 0.001 ***

# m5 - m4 == 0          0.0161790  0.0020367   7.944  < 0.001 ***
# m5 - m1 == 0          0.0226019  0.0020367  11.097  < 0.001 ***
# m5 - m2 == 0          0.0201447  0.0020367   9.891  < 0.001 ***
# poisson - m5 == 0    -0.0191785  0.0020367  -9.416  < 0.001 ***
# m5 - m3 == 0          0.0181480  0.0020367   8.910  < 0.001 ***










# 
# #### PANDIT results
# simulation_topology %>%
#     filter(whichtest == "au") %>%
#     gather(model, pvalue, m1:true) %>%
#     
#     
#     mutate(sig = pvalue < 0.01) %>%
#     filter(sig == TRUE) %>%
#     group_by(name, tree, model) %>%
#     tally() %>%
#     mutate(n = 100* n/20)
#     group_by(model, name, tree, notsig) %>%
#     tally() %>%
#     mutate(perc_in_conf=n/200) %>%
#     mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) -> pandit_au_summary
# 
# 
# 







# 
# 
# pandit_ufb %>% 
#     filter(model != "m1") %>%
#     mutate(model_levels = factor(model, levels=model_levels_nom1, labels = model_labels_nom1),
#         name_levels  = factor(name)) -> ufb_fact
# ufb_fact %>%
#     group_by(model_levels, name_levels) %>%
#     tally() %>%
#     rename(total_considered = n) -> ufb_total_nodes
#      
# ############## False positive nodes ##############
# ufb_fact %>% 
#     mutate(supported = boot >= 95) %>% 
#     filter(supported == TRUE, in_true == FALSE) %>% 
#     group_by(model_levels, name_levels) %>%
#     tally() %>% 
#     complete(model_levels, name_levels, fill = list(n = 0)) %>%
#     left_join(ufb_total_nodes) %>%
#     mutate(percent = n / total_considered) %>%
#     distinct() -> ufb_fp
# 
# ufb_fp %>%
#     ungroup() %>%
#     mutate(controlled = percent <= 0.05) %>%
#     group_by(model_levels) %>%
#     summarize(percent_controlled = 100 * sum(controlled) / 200) -> ufb_fp_controlled
# 
#     
# ufb_fp %>%
#     ggplot(aes(x = model_levels, y = percent, fill = model_levels)) +
#         geom_sina(shape=21, color="grey30", maxwidth=0.75) +
#         #geom_jitter(position=position_jitterdodge(), shape = 21, alpha=0.5,  size=2, color="grey10")  +     
#         scale_fill_manual(values=model_colors_nom1) +         
#         background_grid() + 
#         panel_border() +
#         theme(legend.position = "none",
#               strip.text = element_text(size=8), 
#               axis.text.x = element_text(size=8)) +
#         xlab("Protein Models") + ylab("Percent false positive nodes\nrelative to M1") +
#         geom_hline(yintercept = 0.05, color = "dodgerblue3")  +
#         geom_text(data = ufb_fp_controlled, aes(x = model_levels, y = 0.043, label = paste0(percent_controlled, "%")), fontface="bold", size=3, nudge_x = -0.36) -> pandit_ufb_fp  
# 
# save_plot(paste0(figure_directory,"ufb_pandit_fp.pdf"), pandit_ufb_fp, base_width = 10, base_height=3)
# 
# 
# 
# 
# ##################### Node accuracy #####################################
# give.n <- function(x){
#   return(c(y = median(x) + 0.1, label = length(x)))
# }
# ufb_fact %>% 
#     mutate(supported = boot >= 95) %>% 
#     filter(supported == in_true) %>% 
#     group_by(model_levels, name_levels) %>%
#     tally() %>% 
#     complete(model_levels, name_levels, fill = list(n = 0)) %>%
#     left_join(ufb_total_nodes) %>%
#     mutate(percent = n / total_considered) %>% ## fuck you yy hack
#     distinct() -> ufb_fact_acc
# 
# ufb_fact_acc %>%
#     group_by(model_levels) %>%
#     summarize(meanpercent= 100* round(mean(percent), 3)) -> ufb_fact_acc_means
#  
# ufb_fact_acc %>%   
#     ggplot(aes(x = model_levels, y = percent, fill = model_levels)) +
#         geom_sina(shape=21, color="grey30", maxwidth=0.75) +
#         #geom_jitter(position=position_jitterdodge(), shape = 21, alpha=0.5,  size=2, color="grey10")  +     
#         scale_fill_manual(values=model_colors_nom1) +         
#         background_grid() + 
#         panel_border() +
#         theme(legend.position = "none",
#               strip.text = element_text(size=8), 
#               axis.text.x = element_text(size=8)) +
#         scale_y_continuous(limits=c(0,1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
#         geom_text(data = ufb_fact_acc_means, y=0.9, aes(x = model_levels, label = paste0(meanpercent, "%")), nudge_x = -0.35)+
#         #stat_summary(fun.y = "mean", geom = "text", aes(x = model_levels, label = round(..y.., 2)), position = position_nudge(x = -0.4)) +
#         xlab("Protein Models") + ylab("Percent accurate nodes\nrelative to M1")  -> pandit_ufb_acc 
# 
# save_plot(paste0(figure_directory,"ufb_pandit_accuracy.pdf"), pandit_ufb_acc, base_width = 10, base_height=3)
# 
# 
# ufb_fp %>%
#     filter(name_levels %in% c("HIV", "NP")) %>%
#     ggplot(aes(x = model_levels, y = percent, fill = model_levels)) + 
#         geom_jitter(position=position_jitterdodge(), shape = 21, alpha=0.5,  size=2, color="grey10")  +     
#         scale_fill_manual(values=model_colors, name = "Protein Model") +         
#         facet_grid(name_levels~tree_levels) +
#         background_grid() + 
#         panel_border() +
#         theme(legend.position = "none",
#               strip.text = element_text(size=8), 
#               axis.text.x = element_text(size=6, angle=30)) +
#         xlab("Protein Models") + ylab("Percent false positive nodes") +
#         geom_hline(yintercept = 0.05, color = "dodgerblue3")  -> sim_ufb_fp
# save_plot(paste0(figure_directory,"ufb_simulation_fp_HIV-NP.pdf"), sim_ufb_fp, base_width = 10, base_height=4)
# 
# 
# 
# 
# 

