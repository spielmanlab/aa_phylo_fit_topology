import os
import re
import sys
from copy import deepcopy


def run_save_iqtree(alignment_file, data_type, model, outname, true_tree, threads):
    
    ### Optimize tree under model ###
    if not os.path.exists(outname + "_inferredtree.iqtree"):
        print("Inferring tree", model)
        os.system("iqtree -nt " + str(threads) + " -m " + model + " -s " + alignment_file + " -st " + data_type)
        os.system("mv " + alignment_file + ".log " + outname + "_inferredtree.log")
        os.system("mv " + alignment_file + ".treefile " + outname + "_inferredtree.treefile")
        os.system("mv " + alignment_file + ".iqtree " + outname + "_inferredtree.iqtree")
        os.system("rm " + alignment_file + ".*")    
    
    ### Optimize parameters on true tree ###
    if not os.path.exists(outname + "_optimizedtruetree.iqtree"):
        print("Optimzing true tree", model)
        os.system("iqtree -nt " + str(threads) + " -m " + model + " -s " + alignment_file + " -st " + data_type + " -te " + true_tree)
        os.system("mv " + alignment_file + ".log " + outname + "_optimizedtruetree.log")
        os.system("mv " + alignment_file + ".treefile " + outname + "_optimizedtruetree.treefile")
        os.system("mv " + alignment_file + ".iqtree " + outname + "_optimizedtruetree.iqtree")
        os.system("rm " + alignment_file + ".*")     
      


def main():
    name    = sys.argv[1]
    tree    = sys.argv[2]
    repl    = sys.argv[3]
    threads = sys.argv[4]
    
    rawname = name + "_" + tree + "_rep" + repl + "_"
    
    alignment_path    = "simulations/alignments/"
    selection_path    = "selected_models/"
    fitted_tree_path  = "fitted_trees/" 
    true_tree_path    = "simulations/true_trees/"
    hb_path           = "hb_models/"
    pogomodel_path    = "pogofit_models/"

    with open(selection_path + rawname + "selected_models.csv", "r") as f:
        selection_lines = f.readlines()[1:]

    #print(selection_lines)
    for line in selection_lines:
        splitline = line.strip().split(",")
        # "pandit,tree,repl,datatype,model,kind,ic\n"


        data_type   = splitline[3]
        model       = splitline[4]
        kindofmodel = splitline[5]
        ic          = splitline[6].strip()
    
        
        alignment_file = alignment_path + rawname + data_type + ".fasta"
        true_tree_file = true_tree_path + "tree" + tree + ".tree"
    
        outname_prefix = fitted_tree_path + rawname + data_type + "_"
        print(alignment_file)
        if data_type == "AA" and kindofmodel  == "best": 
            outname = outname_prefix + "poisson"
            run_save_iqtree(alignment_file, data_type, "Poisson", outname, true_tree_file)
    
            hbmodel = hb_path + name + "_HB.paml+G+F"
            outname = outname_prefix + "hbstyle"
            run_save_iqtree(alignment_file, data_type, hbmodel, outname, true_tree_file)

            pogomodel = pogomodel_path + rawname + "AA.dat.POGOFIT.paml+G"  ### already has +F
            outname = outname_prefix + "pogofit"
            run_save_iqtree(alignment_file, data_type, pogomodel, outname, true_tree_file)

        
        
        outname = outname_prefix + "_".join([model,kindofmodel,ic])
        
        ### in case ICs gave the same model, don't reinvent the wheel.
       # run_the_trees = True
        
        #for check_ic in all_ic:
        #    check_outname = outname_prefix + "_".join([model,kindofmodel,check_ic])
        #    if os.path.exists(check_outname + ".log"): 
        #        for suff in [".log", ".treefile"]:
        #            for tree in ["_truetree", "_optimizedtree"]:
        #                os.system("cp " + check_outname + tree + suff + outname + tree + suff)
        #        run_the_trees = False
        #        break   
         
        #if run_the_trees:
        run_save_iqtree(alignment_file, data_type, model, outname, true_tree_file)

main()
