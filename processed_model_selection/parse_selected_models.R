args = commandArgs(trailingOnly=TRUE)
library(tidyverse)

type = args[1]
nreps <- 20
if (type == "empirical"){
    infile <- "processed_model_selection/all_model_selection_empirical.csv"
    outfile <- "processed_model_selection/range_model_selection_empirical.csv"
}
if (type == "pandit"){
    infile <- "processed_model_selection/all_model_selection_pandit.csv"
    outfile <- "processed_model_selection/range_model_selection_pandit.csv"
}
dat <- read_csv(infile)



 

                             

if (type == "pandit")
{
                              
    selected_models <- tibble(name = as.character(),
                          tree = as.character(),
                          repl = as.character(),  ## since NA's
                          bic   = as.numeric(),
                          model = as.character(),
                          modelr = as.integer())
                          
    ## gotta loop since sometimes multiple rows get returned for similarly fitting models
    for (namex in unique(dat$name)){
        dat %>% filter(name == namex) -> subdat
        b <- subdat$bic
        r_step <- (max(b) - min(b)) / 4
        for (r in 1:5) {
            subdat %>%
                # mutate(thisbic = quantile(bic)[q],
#                         diffbic = abs(bic - thisbic)) %>%
#                 filter(diffbic == min(diffbic)) %>%
#                 ungroup() %>%
                mutate(thisbic = min(bic) + (r-1)*r_step, 
                       diffbic = abs(bic - thisbic)) %>%
                filter(diffbic == min(diffbic)) %>%
                dplyr::select(name, model, bic) %>%
                mutate(modelr = r) -> tempdat
                if (nrow(tempdat) >1 )
                {
                    tempdat <- tempdat[1,] %>% mutate("tree" = "NA", "repl" = "NA")
                }
                selected_models <- bind_rows(selected_models, tempdat)
    }}

}

if (type == "empirical")
{

    selected_models <- tibble(name = as.character(),
                          tree = as.character(),
                          repl = as.integer(),  ## since NOT NA's
                          bic   = as.numeric(),
                          model = as.character(),
                          modelr = as.integer())
                          
                           
    ## gotta loop since sometimes multiple rows get returned for similarly fitting models
    for (namex in unique(dat$name)){
        for(treex in unique(dat$tree)){
            for(replx in 1:nreps){
                dat %>% filter(name == namex, tree == treex, repl == replx) -> subdat
                b <- subdat$bic
                r_step <- (max(b) - min(b)) / 4
                for (r in 1:5) {
                    subdat %>%
#                        mutate(thisbic = quantile(bic)[q],
#                               diffbic = abs(bic - thisbic)) %>%
#                        filter(diffbic == min(diffbic)) %>%
 #                       ungroup() %>%
                        mutate(thisbic = min(bic) + (r-1)*r_step, 
                            diffbic = abs(bic - thisbic)) %>%
                        filter(diffbic == min(diffbic)) %>%
                        dplyr::select(name, tree, repl, model) %>%
                        mutate(modelr = r) -> tempdat
                        if (nrow(tempdat) >1 )
                        {
                            tempdat <- tempdat[1,]
                        }
                        selected_models <- bind_rows(selected_models, tempdat)
    }}}}

}
write_csv(selected_models, outfile)        
