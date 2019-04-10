library(tidyverse)

dat <- read_csv("results_sh_rtree.csv")

dat %>% 
    gather(model, pvalue, pogofit:true) %>% 
    filter(tree !="rtree100_bl3_rep1") %>%
    mutate(sig = pvalue <= 0.01) %>% 
    group_by(name, tree, model) %>% 
    tally(sig) %>% 
    filter(n>0) %>% 
    arrange(tree,name) %>%  
    print.data.frame()
#    name                 tree   model  n
# 1  Gal4 rtree100_bl0.75_rep1 poisson  1
# 2  Gal4 rtree100_bl0.75_rep1    true  1
# 3   HIV rtree100_bl0.75_rep1      q5  4
# 4   LAC rtree100_bl0.75_rep1      hb  1
# 5   LAC rtree100_bl0.75_rep1 pogofit  1
# 6   LAC rtree100_bl0.75_rep1      q5  2
# 7   LAC rtree100_bl0.75_rep1    true  1
# 8    NP rtree100_bl0.75_rep1      q5  1
# 9  Gal4  rtree100_bl1.5_rep1    true  1
# 10  HIV  rtree100_bl1.5_rep1 poisson  2
# 11  HIV  rtree100_bl1.5_rep1      q5 16
# 12  HIV  rtree100_bl1.5_rep1    true  1
# 13  LAC  rtree100_bl1.5_rep1 poisson  6
# 14  LAC  rtree100_bl1.5_rep1      q5  1
# 15   NP  rtree100_bl1.5_rep1      q5  3
# 16   NP  rtree100_bl1.5_rep1    true  7


dat %>% 
    gather(model, pvalue, pogofit:true) %>% 
    filter(tree !="rtree100_bl3_rep1") %>%
    mutate(sig = pvalue <= 0.001) %>% 
    group_by(name, tree, model) %>% 
    tally(sig) %>% 
    filter(n>0) %>% 
    arrange(tree,name) %>%  
    print.data.frame()  
#     name                 tree   model  n
# 1  HIV rtree100_bl0.75_rep1      q5  2
# 2  LAC rtree100_bl0.75_rep1      hb  1
# 3  LAC rtree100_bl0.75_rep1 pogofit  1
# 4  LAC rtree100_bl0.75_rep1      q5  1
# 5  LAC rtree100_bl0.75_rep1    true  1
# 6  HIV  rtree100_bl1.5_rep1      q5 14
# 7  LAC  rtree100_bl1.5_rep1 poisson  1
# 8   NP  rtree100_bl1.5_rep1    true  2