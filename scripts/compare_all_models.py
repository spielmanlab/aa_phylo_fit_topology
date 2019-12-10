"""
All m1 models were JTT, HIVb, or WAG. 
Script to compare matrices.
"""
import pyvolve
import numpy as np
from scipy import stats
from pprint import pprint

NORM_I = 0
NORM_J = 1

## These are all the protein models examined in model selection
IQTREE_MODELS = ["JTT", "HIVb", "WAG", "LG", "cpREV", "mtREV", "Dayhoff", "mtMAM", "mtART", "mtZOA", "VT", "rtREV", "DCMut", "PMB", "HIVw", "JTTDCMut", "FLU", "Blosum62", "mtMet", "mtVer", "mtInv"]


def normalize_flatten_matrix(matrix, i, j):
    norm_mat = np.true_divide(matrix, matrix[i][j])
    assert(norm_mat[i][j] == 1.), "Failed normalization"
    return (np.matrix.flatten(norm_mat))
    
    
def main():

    models = {}
    
    outpath = "../results_analysis/csv_files/"
    for model in IQTREE_MODELS:        
        pyvolve_name = model.lower() + "_matrix"            
        if pyvolve_name == "dcmut_matrix":
            pyvolve_name = "dayhoffdcmut_matrix"
        if pyvolve_name == "mtrev_matrix":
            pyvolve_name = "mtrev24_matrix"
    
        matrix = np.array( eval("pyvolve.empirical_matrices." + pyvolve_name) )
        assert(matrix is not None), "bad matrix"
        models[model.strip()] = normalize_flatten_matrix(matrix, NORM_I, NORM_J)

 
    all_models = list(models.keys())
 
    #### First output file: All the flattened rates ####
#     outstring = ",".join(all_models) + "\n"
#     for i in range(len(models[all_models[0]])):
#         line = ""
#         for model in all_models:
#             line += str(models[model][i]) + ","
#         outstring += line.rstrip(",") + "\n"
#     with open(outpath + "all_models_rates.csv", "w") as f:
#         f.write(outstring)

    #### Second output file: matrix correlations ####
    outstring = "model1,model2,r\n"
    done = []

    for m1 in models:
        for m2 in models:
            sorted_str = "".join(sorted([m1, m2]))
            if m1 != m2 and sorted_str not in done:
                done.append (sorted_str)
                r = stats.pearsonr(models[m1], models[m2])[0]
                outstring += ",".join( [m1, m2, str(r)] ) + "\n"
            
    with open(outpath + "all_models_pearson.csv", "w") as f:
        f.write(outstring.strip())        
      
    

 
main()


   
    # 
#     
# jtt_raw  = np.array( pyvolve.empirical_matrices.jtt_matrix )
# hivb_raw = np.array( pyvolve.empirical_matrices.hivb_matrix )
# wag_raw  = np.array( pyvolve.empirical_matrices.wag_matrix )
# 
# 
# models = {"JTT":  normalize_flatten_matrix(jtt_raw, NORM_I, NORM_J), 
#           "HIVb": normalize_flatten_matrix(hivb_raw, NORM_I, NORM_J), 
#           "WAG":  normalize_flatten_matrix(wag_raw, NORM_I, NORM_J)}
# 
# 
# 
# outstring = "JTT,HIVb,WAG\n"
# for i in range(len(models["JTT"])):
#     outstring += ",".join( [str(models["JTT"][i]),
#                            str(models["HIVb"][i]),
#                            str(models["WAG"][i])] )
#     outstring += "\n"
# with open(outpath + "m1_rates.csv", "w") as f:
#     f.write(outstring)
# 
# 
