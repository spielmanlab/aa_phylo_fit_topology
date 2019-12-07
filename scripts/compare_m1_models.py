import pyvolve
import numpy as np
from scipy import stats


NORM_I = 0
NORM_J = 1

outpath = "../results_analysis/"

def normalize_flatten_matrix(matrix, i, j):
    norm_mat = np.true_divide(matrix, matrix[i][j])
    assert(norm_mat[i][j] == 1.), "Failed normalization"
    return (np.matrix.flatten(norm_mat))
    
jtt_raw  = np.array( pyvolve.empirical_matrices.jtt_matrix )
hivb_raw = np.array( pyvolve.empirical_matrices.hivb_matrix )
wag_raw  = np.array( pyvolve.empirical_matrices.wag_matrix )


models = {"JTT":  normalize_flatten_matrix(jtt_raw, NORM_I, NORM_J), 
          "HIVb": normalize_flatten_matrix(hivb_raw, NORM_I, NORM_J), 
          "WAG":  normalize_flatten_matrix(wag_raw, NORM_I, NORM_J)}



outstring = "JTT,HIVb,WAG\n"
for i in range(len(models["JTT"])):
    outstring += ",".join( [str(models["JTT"][i]),
                           str(models["HIVb"][i]),
                           str(models["WAG"][i])] )
    outstring += "\n"
with open(outpath + "m1_rates.csv", "w") as f:
    f.write(outstring)


done = []
outstring = "model1,model2,r,rho,rmse\n"
for m1 in models:
    for m2 in models:
        sorted_str = "".join(sorted([m1, m2]))
        if m1 != m2 and sorted_str not in done:
            done.append (sorted_str)
            r = stats.pearsonr(models[m1], models[m2])[0]
            #rho = stats.spearmanr(models[m1], models[m2])[0]
            rmse = np.sqrt( np.mean( (models[m1] - models[m2])**2 ) )
            outstring += ",".join( [m1, m2, str(r), str(rmse)] ) + "\n"
            
with open(outpath + "m1_comparison.csv", "w") as f:
    f.write(outstring.strip())        
    