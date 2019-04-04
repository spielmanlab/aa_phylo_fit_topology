"""
SJS. Run model selection w/ iqtree on simulations and save csv with best/worst model by each IC for DNA, CODON, and AA models.
"""

import os
import re
import sys


def parse_all_fits(file, prefix):

    start_line = " No. Model         -LnL         df  AIC          AICc         BIC"
    stop_line  = "Akaike Information Criterion"
    with open(file, "r") as f:
        lines = [line for line in f.readlines() if "WARNING" not in line and "Ascertainment bias" not in line]

    for i in range(len(lines)):
        if lines[i].startswith(start_line):
            start_parse = i + 1
        if lines[i].startswith(stop_line):
            stop_parse = i 
            break
    
    #model,logl,df,aic,aicc,bic"
    outstring = ""
    for i in range(start_parse, stop_parse):
        splitline = re.split("\s+", lines[i].strip())
        model = splitline[1]        
        logl  = splitline[2]
        df = splitline[3]
        aic = splitline[4]
        aicc = splitline[5]
        bic = splitline[6]    
        outstring += prefix + "," + ",".join([model, logl, df, aic, aicc, bic]) + "\n"
    return(outstring)    


alignment_path = "simulations/alignments/"
output_path    = "selected_models/" ## mv log files here


outfile = "all_model_selection.csv"
outstring = "name,tree,repl,model,logl,df,aic,aicc,bic\n"
for name in ["NP", "HA", "HIV", "Gal4", "LAC"]:
    for tree in ["rtree100_bl0.3_rep1", "rtree100_bl1.5_rep1", "rtree100_bl3_rep1"]:
        for repli in range(1,21):

            repl = str(repli)
            rawname = name + "_" + tree + "_rep" + repl + "_AA"
            prefix  = ",".join([name, tree, repl])
            print(prefix)
            alignment_file       = alignment_path + rawname + ".fasta"
            model_selection_file = output_path + rawname + ".model_selection_log"
            #print(model_selection_file)
            outstring += parse_all_fits(model_selection_file, prefix)
  
outstring = outstring.strip()
with open(outfile, "w") as f:
    f.write(outstring)

    

