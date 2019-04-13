library(ape)

tnames <- c("anderson", "dosreis", "greenalga", "greenplant", "opisthokonta", "prum", "salichos", "yeast", "ruhfel", "oconnell")

for (t in tnames)
{

    tree <- read.tree(paste0(t, ".tree"))
    tree2 <- multi2di(tree)
    print(t)
    print(length(tree2$tip.label))
    print(sum(tree2$edge.length))
    print("------------------------")
    write.tree(tree2, paste0(t, "_resolved.tree"))
}

# [1] "anderson"
# [1] 179
# [1] 6.617235
# [1] "------------------------"
# [1] "dosreis"
# [1] 274
# [1] 13.87941
# [1] "------------------------"
# [1] "greenalga"
# [1] 23
# [1] 2.481443
# [1] "------------------------"
# [1] "greenplant"
# [1] 360
# [1] 43.05466
# [1] "------------------------"
# [1] "opisthokonta"
# [1] 70
# [1] 20.92042
# [1] "------------------------"
# [1] "prum"
# [1] 200
# [1] 5.211793
# [1] "------------------------"
# [1] "salichos"
# [1] 23
# [1] 9.450554
# [1] "------------------------"
# [1] "yeast"
# [1] 96
# [1] 35.57072
# [1] "------------------------"
# [1] "ruhfel"
# [1] 360
# [1] 24.66852
# [1] "------------------------"
# [1] "oconnell"
# [1] 66
# [1] 4.822812
# [1] "------------------------"
