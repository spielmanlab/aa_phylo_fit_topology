"""
Entropy (in codons) of all simulations per site.
"""
import pyvolve
import numpy as np


def calculate_entropy(f):
    return -1. * np.sum ( f * np.log(f) )

outstring = "name,site,entropy\n"

all_names = ["HA", "NP", "HIV", "1RII", "1IBS", "1R6M"]
dms_sims = ["HA", "NP", "HIV"]
prefpath = "../simulations/preferences/"
for name in all_names:
    print(name)
    prefs = np.loadtxt(prefpath + name + "_prefs.csv", delimiter=",")
    i = 0
    for siteprefs in prefs:
        # preferences are already fitness for yeast
        if name in dms_sims:
            sitefit = np.log(  siteprefs/np.sum(siteprefs)   ) ## renormalize DMS, my tolerance is a bit more stringent
        else:
            sitefit = siteprefs
        m = pyvolve.Model("mutsel", {"fitness": sitefit})
        
        freqs = m.extract_state_freqs()
        h = str(calculate_entropy(freqs))
        outstring += ",".join([name, str(i),h]) + "\n"
        i+=1
with open("../results_analysis/simulation_site_entropy.csv", "w") as f:
    f.write(outstring.strip())