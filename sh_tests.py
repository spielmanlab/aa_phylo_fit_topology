import os
import sys
import re
import pprint
    
type = sys.argv[1]
outfile = "results_sh_" + type + ".csv"


fitted_path = "fitted_trees/"
sim_path    = "simulations/"
true_tree_path = sim_path + "true_trees/"


iqtree_topline = "Tree      logL    deltaL  bp-RELL    p-KH     p-SH       c-ELW\n"

outstring      = "type,name,tree,repl,pogofit,hb,q1,q2,q3,q4,q5,poisson,true\n"
inferencepaths = {"rtree": "fitted_trees/", "empirical": "fitted_trees_empirical/"}
alignmentpaths = {"rtree": sim_path + "alignments/", "empirical": sim_path + "alignments_empirical/"}
quantilefile   = {"rtree": "quantile_model_selection.csv", "empirical": "quantile_model_selection_empirical.csv"}
dms_list       = {"rtree": ["NP", "LAC", "Gal4", "HA", "HIV"], "empirical": ["HA"]}
treenames      = {"rtree": ["rtree100_bl0.3_rep1", "rtree100_bl1.5_rep1", "rtree100_bl0.75_rep1", "rtree100_bl3_rep1"],
                  "empirical": ["anderson", "dosreis", "greenalga", "greenplant", "opisthokonta", "prum", "ruhfel", "salichos", "yeast"]       
                 }
    
reps           = 20


model_order = ["pogofit", "hb", "q1", "q2", "q3", "q4", "q5", "poisson"]



### Get the fit models ###
# name,tree,repl,model,modelq
# NP,rtree100_bl0.3_rep1,1,JTT+F+I+G4,1
print("Determing all best-fitting models")
fitmodels = {}
with open(quantilefile[type], "r") as f:
    qlines = f.readlines()
for qline in qlines:
    name = qline.split(",")[0]
    tree = qline.split(",")[1]
    rep = "rep" + qline.split(",")[2]
    prefix = "_".join([name, tree, rep, "AA"])
    model = qline.split(",")[3]
    quant = qline.split(",")[4].strip()
    if quant == "1":
        fitmodels[prefix] = model



print("Testing galore")
for name in dms_list[type]:
    for tree in treenames[type]:
        if tree != "rtree100_bl3_rep1":
            continue 
        for rep in range(1,reps+1):
            rawname = "_".join([name, tree, "rep" + str(rep), "AA"])
            print(rawname)
            ## Create input file with all trees in order of model_order, and then the true tree
            inferred_trees_raw = [x for x in os.listdir(fitted_path) if x.startswith(rawname) and x.endswith("inferredtree.treefile")]
            inferred_trees_ordered = [] ### order based on replicate, model
            for model in model_order:
                for inftree in inferred_trees_raw:
                    if model in inftree:
                        print(model)
                        with open(fitted_path + inftree, "r") as f:
                            ts = f.read().strip()
                        inferred_trees_ordered.append(ts)
            with open(true_tree_path + tree + ".tree", "r") as f:
                truetree = f.read().strip()
            inferred_trees_ordered.append( truetree )
            with open("treelist.trees", "w") as f:
                for tree in inferred_trees_ordered:
                    f.write(tree +"\n")
    
            ## Call the SH test
            
            alnfile = alignmentpaths[type] + rawname + ".fasta"
            os.system("iqtree -s " + alnfile + " -m " + fitmodels[rawname] + " -z treelist.trees -n 0 -zb 1000 -redo")

            with open(alnfile + ".iqtree", "r") as f:
                output = f.readlines()

            x = 0
            for line in output:
                if line == iqtree_topline:
                    break
                x+=1
    
            start = x + 2
            stop = x + 28
    
            p_sh = []
            for i in range(start, stop):
                p_sh.append( re.split("\s+", output[i])[8] )
            p_sh_string = ",".join(p_sh)

    

            outstring += type + "," + ",".join([name, tree,"rep" + str(rep),p_sh_string]) + "\n"
            print(outstring)
            os.system("rm " + alnfile + ".*")
            assert 1==4

with open(outfile, "w") as f:
    f.write(outstring)    
    
