import os
import sys
import dendropy
from dendropy.calculate import treecompare
import pprint

sim_path       = "../simulations/"
truepath       = sim_path + "true_trees/"
#inferencepath  = "../fitted_trees_simulation/"
inferencepath  = "../fitted_with_ufb/"
dms_list       = ["NP"] #, "LAC", "Gal4", "HA", "HIV"]
treenames      = ["andersen", "dosreis", "opisthokonta", "prum", "ruhfel", "salichos", "rayfinned", "spiralia"]       
reps           = 20


fileinfo_order = ["name", "tree", "rep", "model", "optim"]
fitinfo_order  = ["logl", "k", "AIC", "AICc", "BIC"]

outfile = "inference_results_simulation_ufb.csv"
outstring = ",".join(fileinfo_order) + "," + ",".join(fitinfo_order) +",rf_true,treelength\n"

outfile_boot  = "simulation_ufb_splits.csv"
outstring_boot = "name,tree,rep,model,optim,boot,in_true\n"

iqfiles_all = [x for x in os.listdir(inferencepath) if x.endswith(".iqtree")]

for dms in dms_list: 
    print(dms)
    for repindex in range(1,reps+1): ## files indexed from 1
        print(repindex)
        treespaces = {}
        truetrees = {}
        for treename in treenames:
            ts = dendropy.TaxonNamespace()
            treespaces[treename] = ts
            this_true = dendropy.Tree.get_from_path(truepath + treename +"_resolved.tree", "newick", taxon_namespace = treespaces[treename], rooting='force-unrooted')
            this_true.encode_bipartitions()
            truetrees[treename] = this_true
            
        for this_treename in treenames:
            print("  ",this_treename)
            for optim in ["inferredtree"]:#, "optimizedtruetree"]:    

                prefix = dms + "_" + this_treename.replace(".tree","") + "_rep" + str(repindex) + "_AA"
                iqfiles = [x for x in iqfiles_all if x.startswith(prefix) and optim in x]
                
                for iqfile in iqfiles:
                    name = iqfile.replace(".iqtree","")
                    treefile = name + ".contree"  #".treefile"

                    fileinfo = {"name": None, "treerep": None, "tree": None, "model": None}
                    components  = name.split("_")
                    fileinfo["name"]  = components[0]
                    fileinfo["tree"]  = components[1]
                    fileinfo["rep"]   =  components[2].replace("rep","")
                    fileinfo["model"]   = components[4]
                    fileinfo["optim"] = optim

                    try:
                        inftree = dendropy.Tree.get_from_path(
                                    inferencepath + treefile,
                                    "newick", 
                                    taxon_namespace = treespaces[this_treename]
                                   )
                    except:
                        continue
                    inftree.encode_bipartitions()
                    rf = str( treecompare.symmetric_difference( truetrees[this_treename], inftree) )
                    tl = str(inftree.length())
                    
                    ### bipart ufb compare
                    ### NOTE: UFB are not calculated for POLYTOMIES (a bunch w/ 0 branch lengths)
                    ### These nodes will be IGNORED.
                    bp_prefix = ",".join([fileinfo[x] for x in fileinfo_order])
                    for bp in inftree.bipartition_edge_map: 
                        node = inftree.bipartition_edge_map[bp].head_node 
                        if not node.is_leaf():
                            if node.label: ## near-zero branch lengths aren't included
                                ufb = node.label
                                if bp in truetrees[this_treename].bipartition_edge_map:
                                    outstring_boot += bp_prefix + "," + str(ufb) + ",TRUE\n"
                                else:
                                    outstring_boot += bp_prefix + "," + str(ufb) + ",FALSE\n"
                            else:
                                if node.edge_length is not None:
                                    assert(node.edge_length <= 1e-5), "MISSING BOOTSTRAP?!"

                    #assert 1==3

            
                    fitinfo = {"logl": None, "k": None, "AIC": None, "AICc": None, "BIC": None}
                    with open(inferencepath + iqfile, "r") as f:
                        iqlines = f.readlines()
                    for line in iqlines:
                        if line.startswith("Log-likelihood of the tree"):
                            fitinfo["logl"] = line.split(":")[1].split(" ")[1]
                        if line.startswith("Number of free parameters"):
                            fitinfo["k"] = line.split(" ")[-1].strip()
                        if line.startswith("Akaike information criterion"):
                            fitinfo["AIC"] = line.split(" ")[-1].strip()
                        if line.startswith("Corrected Akaike information criterion "):
                            fitinfo["AICc"] = line.split(" ")[-1].strip()
                        if line.startswith("Bayesian information criterion"):
                            fitinfo["BIC"] = line.split(" ")[-1].strip()
                    part = ",".join([fileinfo[x] for x in fileinfo_order]) + "," + ",".join([fitinfo[x] for x in fitinfo_order]) + "," + rf + "," + tl  + "\n"
                    outstring += part

with open(outfile, "w") as f:
    f.write(outstring.strip())
    
with open(outfile_boot, "w") as f:
    f.write(outstring_boot.strip())
    
    
