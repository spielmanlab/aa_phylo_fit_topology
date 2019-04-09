args = commandArgs(trailingOnly=TRUE)
library(tidyverse)

type = args[1]
nreps <- 20
if (type == "empirical"){
    infile <- "all_model_selection_empirical.csv"
    outfile <- "quantile_model_selection_empirical.csv"
}
if (type == "rtree"){
    infile <- "all_model_selection.csv"
    outfile <- "quantile_model_selection.csv"
}
dat <- read_csv(infile)

selected_models <- tibble(name = as.character(),
                          tree = as.character(),
                          repl = as.integer(),
                          model = as.character(),
                          modelq = as.integer())




## gotta loop since sometimes multiple rows get returned for similarly fitting models
for (namex in unique(dat$name)){
    for(treex in unique(dat$tree)){
        for(replx in 1:nreps){
            dat %>% filter(name == namex, tree == treex, repl == replx) -> subdat
            for (q in 1:5) {
                subdat %>%
                    mutate(thisbic = quantile(bic)[q],
                           diffbic = abs(bic - thisbic)) %>%
                    filter(diffbic == min(diffbic)) %>%
                    ungroup() %>%
                    dplyr::select(name, tree, repl, model) %>%
                    mutate(modelq = q) -> tempdat
                    if (nrow(tempdat) >1 )
                    {
                        tempdat <- tempdat[1,]
                    }
                    selected_models <- bind_rows(selected_models, tempdat)
}}}}
write_csv(selected_models, outfile)        
