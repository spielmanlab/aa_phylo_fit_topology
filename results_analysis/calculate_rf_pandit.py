import os
import sys
import dendropy
from dendropy.calculate import treecompare
import pprint
truepath       = "../pandit_aa_alignments/"
inferencepath = "../fitted_trees_pandit/"

fileinfo_order = ["name", "model", "optim"]
fitinfo_order  = ["logl", "k", "AIC", "AICc", "BIC"]

outstring = ",".join(fileinfo_order) + "," + ",".join(fitinfo_order) +",rf_q1,rfw_q1,treelength\n"


pandits  = [x.replace(".fasta", "") for x in os.listdir(truepath) if x.endswith("fasta")]

for pandit in pandits:
    print(pandit)
    ts = dendropy.TaxonNamespace()
    q1_tree   = dendropy.Tree.get_from_path(inferencepath + pandit + "_q1_inferredtree.treefile" , "newick", taxon_namespace = ts, rooting='force-unrooted')
    
    iqfiles = [x for x in os.listdir(inferencepath) if x.endswith(".iqtree") and x.startswith(pandit)]

   # print(pandit, len(iqfiles))
    for iqfile in iqfiles:
        fileinfo = iqfile.split("_")
        model = fileinfo[1]
        optim = fileinfo[2].split(".")[0]
        treefile = iqfile.replace(".iqtree", ".treefile")
        inftree = dendropy.Tree.get_from_path(
                    inferencepath + treefile,
                    "newick", 
                    taxon_namespace = ts)
        rf_q1 = str( treecompare.symmetric_difference( q1_tree, inftree) )
        rfw_q1 = str( treecompare.weighted_robinson_foulds_distance( q1_tree, inftree) )
        tl = str(inftree.length())

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
        part = ",".join([pandit, model, optim]) + "," + ",".join([fitinfo[x] for x in fitinfo_order]) + "," + rf_q1 + "," + rfw_q1 + "," + tl  + "\n"
        outstring += part
     
with open("inference_results_pandit.csv", "w") as f:
    f.write(outstring.strip())
    
    
