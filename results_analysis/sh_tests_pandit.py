import os
import sys
import re
import pprint
    
    
alignmentpath = "../pandit_aa_alignments/"
inferencepath = "../fitted_trees_pandit/"
quantilefile   = "../processed_model_selection/quantile_model_selection_pandit.csv"

iqtree_topline = ['Tree', 'logL', 'deltaL', 'bp-RELL', 'p-KH', 'p-SH', 'c-ELW']
model_order    = ["q1", "q2", "q3", "q4", "q5", "poisson"]

outfile   = "results_sh_pandit.csv"
outstring = "name,q1,q2,q3,q4,q5,poisson,pandit\n"



### Get the fit models ###
# name,tree,repl,model,modelq
# NP,rtree100_bl0.3_rep1,1,JTT+F+I+G4,1
print("Determing all best-fitting models")
fitmodels = {}
with open(quantilefile, "r") as f:
    qlines = f.readlines()
for qline in qlines:
    name = qline.split(",")[0]
    model = qline.split(",")[3]
    quant = qline.split(",")[4].strip()
    if quant == "1":
        fitmodels[name] = model



print("Testing galore")
for name in list(fitmodels.keys()):
    print(name)
    inferred_trees_raw = [x for x in os.listdir(inferencepath) if x.startswith(name) and x.endswith("inferredtree.treefile")]
    inferred_trees_ordered = [] 
    for model in model_order:
        for inftree in inferred_trees_raw:
            if model in inftree:
                with open(inferencepath + inftree, "r") as f:
                    ts = f.read().strip()
                inferred_trees_ordered.append(ts)
    
    with open(alignmentpath + name + ".tree", "r") as f:
       pandittree = f.read().strip()
    inferred_trees_ordered.append( pandittree )
        
    with open("treelist.trees", "w") as f:
        for ts in inferred_trees_ordered:
           f.write(ts +"\n")
    ## Call the SH test
    alnfile = alignmentpath + name + ".fasta"
    os.system("iqtree -quiet -nt 4 -s " + alnfile + " -m " + fitmodels[name] + " -z treelist.trees -n 0 -zb 1000 -redo")

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
        #print(output[i])
        p_sh.append( re.split("\s+", output[i])[8] )
    p_sh_string = ",".join(p_sh)
    #print(p_sh_string)
        
    outstring += ",".join([name,p_sh_string]) + "\n"
    os.system("rm " + alnfile + ".*")

with open(outfile, "w") as f:
    f.write(outstring)    
    
