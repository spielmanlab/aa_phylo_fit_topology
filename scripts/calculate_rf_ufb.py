import os
import sys
import dendropy
from dendropy.calculate import treecompare
import pprint

THRESHOLD=95 ## calling positive, negative



def obtain_fit_info(iqfile, fitinfo_order):
    fitinfo = {"logl": None, "k": None, "AIC": None, "AICc": None, "BIC": None}
    with open(iqfile, "r") as f:
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
    return ",".join([fitinfo[x] for x in fitinfo_order])
 
 
def compare_ufb(bp_prefix, inftree, truetree):
    ### NOTE: UFB are not calculated for POLYTOMIES
    ### These nodes will be IGNORED.
    inftree.encode_bipartitions()
    truetree.encode_bipartitions()
    ufb_string = ""
    for bp in inftree.bipartition_edge_map: 
        node = inftree.bipartition_edge_map[bp].head_node 
        if not node.is_leaf():
            if node.label: ## near-zero branch lengths aren't included
                ufb = node.label
                ufb_float = float(ufb)
                # TRUE SPLIT
                if bp in truetree.bipartition_edge_map:
                    ufb_string += bp_prefix + "," + str(ufb)  + ",TRUE"
                    if ufb_float >= THRESHOLD:  # true split and high boot
                        ufb_string += ",TP\n"
                    elif ufb_float < THRESHOLD:  # true split and low boot
                        ufb_string += ",FN\n"
                # FALSE SPLIT
                else:
                    ufb_string += bp_prefix + "," + str(ufb) + ",FALSE"
                    if ufb_float >= THRESHOLD:  # false split and high boot
                        ufb_string += ",FP\n"
                    elif ufb_float < THRESHOLD:  # false split and low boot
                        ufb_string += ",TN\n"
            else:
                if node.edge_length is not None:
                    assert(node.edge_length <= 1e-5), "MISSING BOOTSTRAP?!"
    return ufb_string





type = sys.argv[1]
outpath = "../results_analysis/"
fitinfo_order    = ["logl", "k", "AIC", "AICc", "BIC"]



if type ==  "pandit":

    alignmentpath = "../pandit_aa_alignments/"
    inferencepath = "../fitted_trees_ufb_pandit/"
    fileinfo_order = ["name", "model"]
    outfile_rf  = outpath + "rf_" + type + ".csv"
    outfile_fit = outpath + "fit_" + type + ".csv"
    outstring_fit = "name,model," + ",".join(fitinfo_order) + ",tl\n"
    outstring_rf = "name,model1,model2,rf\n"
    model_order = ["m1", "m2", "m3", "m4", "m5", "poisson", "GTR20"]

    pandits  = [x.replace(".fasta", "") for x in os.listdir(alignmentpath) if x.endswith("fasta")]
    for pandit in pandits:
        print(pandit)
        ts = dendropy.TaxonNamespace()
        trees = []
        
        # collect the trees
        for model in model_order:
            treefile = pandit + "_" + model + "_inferredtree.contree"
            inftree = dendropy.Tree.get_from_path(
                                        inferencepath + treefile,
                                        "newick", 
                                        taxon_namespace = ts, 
                                        rooting='force-unrooted')
            inftree.encode_bipartitions()
            trees.append(inftree)
        # rf and fits
        for i in range(len(model_order)):
            this_model = model_order[i]
            this_tree = trees[i]
            for j in range(i + 1, len(model_order)):
                rf = str(treecompare.symmetric_difference(this_tree, trees[j]))
                outstring_rf += pandit + "," + this_model + "," + model_order[j] + "," + rf + "\n"
        
            iqfile = inferencepath + pandit + "_" + this_model + "_inferredtree.iqtree"
            tl = str(trees[i].length())
            outstring_fit += pandit + "," + this_model + "," + obtain_fit_info(iqfile,fitinfo_order) + "," + tl + "\n"


    with open(outfile_fit, "w") as f:
        f.write(outstring_fit.strip())

    with open(outfile_rf, "w") as f:
        f.write(outstring_rf.strip())        



if type == "simulation":
    alignmentpath  = "../simulations/alignments/"
    inferencepath  = "../fitted_trees_ufb_simulation/" 
    true_tree_path = "../simulations/true_trees/"
    sim_list       = ["LAC", "NP", "HA", "HIV"]
    treenames      = ["dosreis","andersen", "opisthokonta", "prum", "ruhfel", "salichos", "rayfinned", "spiralia"]        
    reps           = 20
    fileinfo_order = ["name", "tree", "rep", "model"]
    outstring_rf_fit = ",".join(fileinfo_order) + "," + ",".join(fitinfo_order) +",rf,treelength\n"
    outstring_boot   = ",".join(fileinfo_order) + ",boot,in_true,classif\n"
    outfile_rf_fit  = outpath + "rf_fit_" + type + ".csv"
    outfile_boot   = outpath + "ufb_splits_" + type + ".csv"


    iqfiles_all = [x for x in os.listdir(inferencepath) if x.endswith(".iqtree")]

    for sim in sim_list: 
        print(sim)
        for repindex in range(1,reps+1): ## files indexed from 1
            print(repindex)
            treespaces = {}
            truetrees = {}
            for treename in treenames:
                ts = dendropy.TaxonNamespace()
                treespaces[treename] = ts
                this_true = dendropy.Tree.get_from_path(true_tree_path + treename +"_resolved.tree", "newick", taxon_namespace = treespaces[treename], rooting='force-unrooted')
                this_true.encode_bipartitions()
                                
                truetrees[treename] = this_true
            
            for this_treename in treenames:
                print("  ",this_treename)

                prefix = sim + "_" + this_treename.replace(".tree","") + "_rep" + str(repindex) + "_AA"
                iqfiles = [x for x in iqfiles_all if x.startswith(prefix)]
                assert(len(iqfiles) == 7), "bad iqfiles"
                for iqfile in iqfiles:
                    name = iqfile.replace(".iqtree","")
                    treefile = name + ".contree"

                    fileinfo = {"name": None, "treerep": None, "tree": None, "model": None}
                    components  = name.split("_")
                    fileinfo["name"]  = components[0]
                    fileinfo["tree"]  = components[1]
                    fileinfo["rep"]   =  components[2].replace("rep","")
                    fileinfo["model"]   = components[4]

                    inftree = dendropy.Tree.get_from_path(
                                    inferencepath + treefile,
                                    "newick", 
                                    taxon_namespace = treespaces[this_treename]
                                   )
                    inftree.encode_bipartitions()

                    rf = str( treecompare.symmetric_difference( truetrees[this_treename], inftree) )
                    tl = str(inftree.length())
                

                    bp_prefix = ",".join([fileinfo[x] for x in fileinfo_order])
                    outstring_boot += compare_ufb(bp_prefix, inftree, truetrees[this_treename])

                    fit_string = obtain_fit_info(inferencepath + iqfile, fitinfo_order)

                    part = ",".join([fileinfo[x] for x in fileinfo_order]) + "," + fit_string + "," + rf + "," + tl  + "\n"
                    outstring_rf_fit += part

    with open(outfile_rf_fit, "w") as f:
        f.write(outstring_rf_fit.strip())

    with open(outfile_boot, "w") as f:
        f.write(outstring_boot.strip())








