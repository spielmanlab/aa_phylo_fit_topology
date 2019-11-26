"""
Entropy (in codons) of all simulations per site.
"""
import pyvolve
import numpy as np
g = pyvolve.Genetics()


def codon_freqs_to_aa_freqs(codonfreqs):
    '''
        Codon codon frequencies to amino-acid frequencies.
    '''
    cf = dict(zip(g.codons, codonfreqs))
    aa_freqs = []
    for family in g.genetic_code:
        total = 0.
        for syn in family:
            total += cf[syn]
        aa_freqs.append(total)
    aa_freqs = np.array(aa_freqs)
    assert((1. - np.sum(aa_freqs)) <= 1e-8), "\nImproperly converted codon to amino acid frequencies."
    return aa_freqs



def calculate_entropy(f):
    return -1. * np.sum ( f * np.log(f) )

outstring = "name,site,entropy\n"

all_names = ["LAC", "NP", "HA", "HIV"]
prefpath = "../simulations/preferences/"
for name in all_names:
    print(name)
    prefs = np.loadtxt(prefpath + name + "_prefs.csv", delimiter=",")
    i = 0
    for siteprefs in prefs:
        sitefit = np.log(  siteprefs/np.sum(siteprefs)   ) ## renormalize DMS, my tolerance is a bit more stringent

        m = pyvolve.Model("mutsel", {"fitness": sitefit})
        
        codon_freqs = m.extract_state_freqs()
        freqs    = codon_freqs_to_aa_freqs(codon_freqs)
        h = str(calculate_entropy(freqs))
        outstring += ",".join([name, str(i),h]) + "\n"
        i+=1
with open("../results_analysis/simulation_site_entropy.csv", "w") as f:
    f.write(outstring.strip())