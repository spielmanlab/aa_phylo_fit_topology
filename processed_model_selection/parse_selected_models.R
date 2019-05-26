args = commandArgs(trailingOnly=TRUE)
library(tidyverse)

type = args[1]
nreps <- 20
if (type == "simulation"){
    infile <- "processed_model_selection/all_model_selection_simulation.csv"
    outfile <- "processed_model_selection/quantile_model_selection_simulation.csv"
}
if (type == "pandit"){
    infile <- "processed_model_selection/all_model_selection_pandit.csv"
    outfile <- "processed_model_selection/quantile_model_selection_pandit.csv"
}
dat <- read_csv(infile)


if (type == "pandit")
{
                              
    selected_models <- tibble(name = as.character(),
                             bic   = as.numeric(),
                             model_name = as.character(),
                             modelm = as.integer())
                          
    ## gotta loop since sometimes multiple rows get returned for similarly fitting models
    for (namex in unique(dat$name)){
        dat %>% filter(name == namex) -> subdat
        b <- subdat$bic
        for (m in 1:5) {
            subdat %>%
                mutate(thisbic = quantile(bic)[m],
                        diffbic = abs(bic - thisbic)) %>%
                filter(diffbic == min(diffbic)) %>%
                ungroup() %>%
                dplyr::select(name, model, bic) %>%
                rename(model_name = model) %>%
                mutate(modelm = m) -> tempdat
                if (nrow(tempdat) >1 ) tempdat <- tempdat[1,]
                selected_models <- bind_rows(selected_models, tempdat)
    }}

}

if (type == "simulation")
{

    selected_models <- tibble(name = as.character(),
                          tree = as.character(),
                          repl = as.integer(),  ## since NOT NA's
                          bic   = as.numeric(),
                          model_name = as.character(),
                          modelm = as.integer())
                          
                           
    ## gotta loop since sometimes multiple rows get returned for similarly fitting models
    for (namex in unique(dat$name)){
        for(treex in unique(dat$tree)){
            for(replx in 1:nreps){
                dat %>% filter(name == namex, tree == treex, repl == replx) -> subdat
                b <- subdat$bic
                for (m in 1:5) {
                    subdat %>%
                        mutate(thisbic = quantile(bic)[m],
                               diffbic = abs(bic - thisbic)) %>%
                        filter(diffbic == min(diffbic)) %>%
                        ungroup() %>%
                        dplyr::select(name, tree, repl, model, bic) %>%
                        rename(model_name = model) %>%
                        mutate(modelm = m) -> tempdat
                        if (nrow(tempdat) >1 ) tempdat <- tempdat[1,]
                        selected_models <- bind_rows(selected_models, tempdat)
    }}}}

}
write_csv(selected_models, outfile)        
