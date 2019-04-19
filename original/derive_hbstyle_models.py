import os
import sys
import pyvolve
import numpy as np
 
"""
    This script calculates an ANALOGY for the protein exchangability model.
    
    Exchangeability model can be decomposed as Qij = Rij*F_j, where Rij is an exchangeability ("rate") and F_j is the target frequency.
    For a given alignment, swMutSel has inferred site-wise fitness and global mutation rates. We can use these to get site-wise Rij. 
    For making the models we will put 0.05's all around and just do +F for inference.   
    
    Save Q matrix (with 0 diagonal) in text format.
    Save PAML model (decomposed rates + dummy 0.05 freqs from the HB model) to text.
"""


ZERO    = 1e-8
G            = pyvolve.Genetics()
GENETIC_CODE = G.genetic_code
CODONS       = G.codons
PREFPATH   = "simulations/preferences/"
OUTPATH    = "hb_models/"
PAML_AA    = "ARNDCQEGHILKMFPSTWYV"
PYVOLVE_AA = "ACDEFGHIKLMNPQRSTVWY"


def codon_freqs_to_aa_freqs(codonfreqs):
    '''
        Codon codon frequencies to amino-acid frequencies.
    '''
    g = pyvolve.Genetics()
    cf = dict(zip(g.codons, codonfreqs))
    aa_freqs = []
    for family in g.genetic_code:
        total = 0.
        for syn in family:
            total += cf[syn]
        aa_freqs.append(total)
    aa_freqs = np.array(aa_freqs)
    assert((1. - np.sum(aa_freqs)) <= ZERO), "\nImproperly converted codon to amino acid frequencies."
    return aa_freqs

mu = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}
paml_frequencies = "\n" +  " ".join(["0.05"]*20)

def main():    

    pref_files = [x for x in os.listdir(PREFPATH) if x.endswith("_prefs.csv")]
    for file in pref_files:
        name = file.split("_")[0]
        outfilepaml = OUTPATH + name + "_HB.paml" 
        prefs = np.loadtxt(PREFPATH + file, delimiter=",")

        mean_codon_matrix = np.zeros([61,61])
        for siteprefs in prefs:
            sitefit = np.log(  siteprefs/np.sum(siteprefs)   ) ## renormalize, my tolerance is a bit more stringent
            m = pyvolve.Model("mutsel", {"fitness": sitefit, "mu": mu})
            mean_codon_matrix += m.extract_rate_matrix() ## 61x61
        mean_codon_matrix /= len(prefs)
        
        
        ## Convert to amino acid model assuming equal frequencies
        aa_matrix = np.zeros([20,20])
        for i in range(20):
            source_codons = GENETIC_CODE[i]
            for j in range(i+1, 20):
                target_codons = GENETIC_CODE[j]
                
                aa_rate = 0.
                x = 0
                for s in source_codons:
                    source_index = CODONS.index(s)
                    for t in target_codons:
                        target_index = CODONS.index(t)
                        this_rate = mean_codon_matrix[source_index][target_index]
                        assert(abs(this_rate) == this_rate), "rate is negative somehow FREAK OUT" 
                        aa_rate += this_rate
                        x += 1
                aa_rate /= x
                
                aa_matrix[i][j] = aa_rate
                aa_matrix[j][i] = aa_rate
        for i in range(20):
            aa_matrix[i][i] = -1. * np.sum(aa_matrix[i])  
            
        ## Decompose to get rate matrix (not q matrix)
        decomped_aa_matrix = np.dot(aa_matrix, np.diag(1.0 / np.repeat(0.05, 20)))  ## frequencies will be equal based on how calculations above were done
        assert(np.allclose(decomped_aa_matrix, decomped_aa_matrix.T)), "\n decomp'd rate matrix not symmetric"

        ## Save to PAML format
        paml = ""
        full = np.zeros([20,20])
        for s in range(20):
            for t in range(20):
                if s==t:
                    continue
                saa = PAML_AA[s] ## paml
                taa = PAML_AA[t] ## paml
                full[s][t] = decomped_aa_matrix[ PYVOLVE_AA.index(saa) ][ PYVOLVE_AA.index(taa) ]
        assert(np.allclose(full, full.T)), "\n paml converting rate matrix not symmetric"
        for i in range(20):
            x = 0
            row = ""
            while x < i:
                row += str(full[x][i]) + " " 
                x+=1
            paml += row.strip() + "\n"

        paml += paml_frequencies
        with open(outfilepaml, "w") as f:
            f.write(paml.strip())

main()

   
