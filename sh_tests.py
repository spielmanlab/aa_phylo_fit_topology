import os
import sys
import re
import pprint
    


sim_path    = "simulations/"
true_tree_path = sim_path + "true_trees/"

inferencepaths = "fitted_trees_empirical/
alignmentpath = sim_path + "alignments_empirical/
quantilefile   = "quantile_model_selection_empirical.csv"
dms_list       = ["NP", "LAC", "Gal4", "HA", "HIV"]
treenames      = ["anderson", "dosreis", "greenalga", "opisthokonta", "prum", "ruhfel", "salichos", "rayfinned", "spiralia"]       
reps           = 20


iqtree_topline = ['Tree', 'logL', 'deltaL', 'bp-RELL', 'p-KH', 'p-SH', 'c-ELW']

model_order    = ["pogofit", "hb", "q1", "q2", "q3", "q4", "q5", "poisson"]

outfile   = "results_sh_empirical.csv"
outstring = "name,tree,repl,q1,q2,q3,q4,q5,poisson,true\n"



### Get the fit models ###
# name,tree,repl,model,modelq
# NP,rtree100_bl0.3_rep1,1,JTT+F+I+G4,1
print("Determing all best-fitting models")
fitmodels = {}
with open(quantilefile, "r") as f:
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
for name in dms_list:
    print(name)
    for tree in treenames:
        print("  ", tree)
        for rep in range(1,reps+1):
            print("  ", rep)
            rawname = "_".join([name, tree, "rep" + str(rep), "AA"])
            #print(rawname)
            ## Create input file with all trees in order of model_order, and then the true tree
            inferred_trees_raw = [x for x in os.listdir(inferencepath) if x.startswith(rawname) and x.endswith("inferredtree.treefile")]
            inferred_trees_ordered = [] ### order based on replicate, model
            for model in model_order:
                for inftree in inferred_trees_raw:
                    if model in inftree:
                        #print(model)
                        with open(inferencepath + inftree, "r") as f:
                            ts = f.read().strip()
                        inferred_trees_ordered.append(ts)
            with open(true_tree_path + tree + ".tree", "r") as f:
                truetree = f.read().strip()
            inferred_trees_ordered.append( truetree )
            with open("treelist.trees", "w") as f:
                for ts in inferred_trees_ordered:
                    f.write(ts +"\n")
    
            ## Call the SH test
            alnfile = alignmentpath + rawname + ".fasta"
            os.system("iqtree -quiet -nt 4 -s " + alnfile + " -m " + fitmodels[rawname] + " -z treelist.trees -n 0 -zb 1000 -redo")

            ## Parse out p-values
            with open(alnfile + ".iqtree", "r") as f:
                output = f.readlines()
            x = 0
            for line in output:
                if re.split("\s+", line.strip()) == iqtree_topline:
                    break
                x+=1
            start = x + 2
            stop = x + len(model_order) + 3
            p_sh = []
            for i in range(start, stop):
                p_sh.append( re.split("\s+", output[i])[8] )
            p_sh_string = ",".join(p_sh)


            outstring += ",".join([name, tree, str(rep),p_sh_string]) + "\n"
            #print(outstring)
            os.system("rm " + alnfile + ".*")

with open(outfile, "w") as f:
    f.write(outstring)    
    
