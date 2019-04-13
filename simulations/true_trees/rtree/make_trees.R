library(ape)


ntaxa <- 100

for (i in 1:10){
    for (bl in c(0.3, 1.5, 3.0)){
        
        t <- rtree(ntaxa)
        t <- compute.brlen(t, bl)
        name <- paste0("rtree", ntaxa, "_bl",bl, "_rep",i, ".tree")
        write.tree(t, name)
    }
}