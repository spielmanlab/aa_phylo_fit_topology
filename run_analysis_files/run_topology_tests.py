import os
import sys
import re



iqtree_topline = ['Tree', 'logL', 'deltaL', 'bp-RELL', 'p-KH', 'p-SH', 'c-ELW']
sh_index = iqtree_topline.index("p-SH")
model_order    = ["m1", "m2", "m3", "m4", "m5", "poisson", "GTR20"]

type = sys.argv[1]
threads = sys.argv[2]
outfilesh        = "shtests_" + type + ".csv"
outfileau        = "autests_" + type + ".csv"
treelist_name  = type + ".trees"
shpath        = "../topology_tests_" + type + "/"
quantilefile  = "../processed_model_selection/quantile_model_selection_" + type + ".csv"


if type == "pandit":
    alignmentpath = "../pandit_aa_alignments/"
    inferencepath = "../fitted_trees_pandit/"
elif type == "simulation":
    alignmentpath  = "../simulations/alignments/"
    inferencepath  = "../fitted_with_ufb/" #"../fitted_trees_simulation/"
    true_tree_path = "../simulations/true_trees/"
    dms_list       = ["NP", "HA", "HIV"]
    treenames      = ["andersen", "dosreis", "opisthokonta", "prum", "ruhfel", "salichos", "rayfinned", "spiralia"]        
    reps           = 20


def determine_best_models(type, quantilefile):
    fitmodels = {}
    with open(quantilefile, "r") as f:
        qlines = f.readlines()
    
    if type == "pandit":
        for qline in qlines:
            name = qline.split(",")[0]
            model = qline.split(",")[2]
            quant = qline.split(",")[3].strip()
            if quant == "1":
                fitmodels[name] = model
    if type == "simulation":
        for qline in qlines:
            name = qline.split(",")[0]
            tree = qline.split(",")[1]
            rep = "rep" + qline.split(",")[2]
            prefix = "_".join([name, tree, rep, "AA"])
            model = qline.split(",")[4]
            quant = qline.split(",")[5].strip()
            if quant == "1":
                fitmodels[prefix] = model
    return(fitmodels)


def parse_tests(output, name):

    x = 0
    for line in output:
        if re.split("\s+", line.strip()) == iqtree_topline:
            break
        x+=1
    start = x + 2
    if type == "simulation":
        stop = x + len(model_order) + 3
    if type == "pandit":
        stop = x + len(model_order) + 2
    p_sh = []
    for i in range(start, stop):
        #print(i, output[i])
        clean = re.split("\s+", output[i].replace("-","").replace("+","").strip())
        p_sh.append( clean[sh_index] )
    p_sh_string = ",".join(p_sh)

    outstring = ",".join([name, p_sh_string]) + "\n"
    
    return(outstring)



def loop_over_tests(type):

    if type == "pandit":
        outstring = "name,m1,m2,m3,m4,m5,poisson,GTR20\n"
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
                with open(treelist_name, "w") as f:
                    for ts in inferred_trees_ordered:
                        f.write(ts +"\n")
            m1tree = inferencepath + name + "_m1_inferredtree.treefile"
            alnfile = alignmentpath + name + ".fasta"
            
            ## Call the SH test
            sh_outfile = shpath + name + ".topology_tests"
            if not os.path.exists(sh_outfile):
                os.system("iqtree -quiet -nt " + threads + " -s " + alnfile + " -m " + fitmodels[name] + " -z " + treelist_name + " -zb 10000 -au -te " + m1tree + " -redo")
                os.system("mv " + alnfile + ".iqtree " + sh_outfile)
                os.system("rm " + alnfile + ".*")
 #           with open(sh_outfile, "r") as f:
 #               lines = f.readlines()
 #           this_out = parse_tests(lines, name)
 #           outstring += this_out
 #           with open(name + ".shtest", "w") as f:
 #               f.write(this_out)


    if type == "simulation":
        outstring  = "name,tree,repl,m1,m2,m3,m4,m5,poisson,GTR20,true\n"
        for name in dms_list:
            print(name)
            for tree in treenames:
                print("  ", tree)
                for rep in range(1,21):
                    print("  ", rep)
                    rawname = "_".join([name, tree, "rep" + str(rep), "AA"])
                    inferred_trees_raw = [x for x in os.listdir(inferencepath) if x.startswith(rawname) and x.endswith("inferredtree.treefile")]
                    inferred_trees_ordered = [] ### order based on replicate, model
                    for model in model_order:
                        for inftree in inferred_trees_raw:
                            if model in inftree:
                                with open(inferencepath + inftree, "r") as f:
                                    ts = f.read().strip()
                                inferred_trees_ordered.append(ts)
                    
                    with open(true_tree_path + tree + "_resolved.tree", "r") as f:
                        truetree = f.read().strip()
                    inferred_trees_ordered.append( truetree )
                    
                    with open(treelist_name, "w") as f:
                        for ts in inferred_trees_ordered:
                            f.write(ts +"\n")
                    
                    alnfile = alignmentpath + rawname + ".fasta"
                    m1tree = inferencepath + rawname + "_m1_inferredtree.treefile"
                    
                    ## Call the SH test
                    sh_outfile = shpath + rawname + ".topology_tests"
                    if not os.path.exists(sh_outfile):
                        os.system("iqtree -quiet -nt " + threads + " -s " + alnfile + " -m " + fitmodels[rawname] + " -z " + treelist_name + " -zb 10000 -au -te " + m1tree + " -redo")
                        os.system("mv " + alnfile + ".iqtree " + sh_outfile)
                        os.system("rm " + alnfile + ".*")
                   # with open(sh_outfile, "r") as f:
                   #     lines = f.readlines()
                   # this_out = parse_tests(lines, name + "," + tree + "," + str(rep))
                   # outstring += this_out
                   # with open(rawname + ".shtest", "w") as f:
                   #     f.write(this_out)
    #return(outstring)

            
        
   
print("Determing all best-fitting models")
fitmodels = determine_best_models(type, quantilefile)
    
print("Run the tests")  
outstring = loop_over_tests(type)
#with open(outfile, "w") as f:
#    f.write(outstring)    
#os.system("rm " + treelist_name)

