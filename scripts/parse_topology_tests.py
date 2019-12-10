import os
import sys
import re

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
    p_au = []
    for i in range(start, stop):
        #print(i, output[i])
        clean = re.split("\s+", output[i].replace("-","").replace("+","").strip())
        p_sh.append( clean[sh_index] )
        p_au.append( clean[au_index] )
    p_sh_string = "sh," + ",".join(p_sh)
    p_au_string = "au," + ",".join(p_au)

    outstring = ",".join([name, p_sh_string]) + "\n" + ",".join([name, p_au_string]) + "\n"
    
    return(outstring)


# 
# Tree      logL    deltaL  bp-RELL    p-KH     p-SH       c-ELW       p-AU
# -------------------------------------------------------------------------
#   1 -32469.16223       0  0.0658 +  0.518 +      1 +     0.121 +    0.757 + 
#   2  -32469.7164 0.55417   0.119 +  0.358 +  0.874 +     0.114 +    0.471 + 
#   3 -32469.44438 0.28214   0.321 +  0.443 +  0.816 +     0.272 +    0.644 + 
#   4 -32471.97116  2.8089   0.131 +  0.234 +  0.622 +     0.106 +     0.29 + 
#   5 -32485.71737  16.555  0.0678 + 0.0989 +  0.134 +    0.0628 +   0.0963 + 
#   6 -32469.58521 0.42298  0.0939 +  0.347 +  0.907 +    0.0996 +    0.449 + 
#   7 -32469.16277 0.00053711  0.0963 +  0.481 +  0.974 +     0.121 +    0.755 + 
#   8 -32469.44838 0.28615   0.105 +  0.315 +   0.92 +     0.103 +    0.422 + 

iqtree_topline = ['Tree', 'logL', 'deltaL', 'bp-RELL', 'p-KH', 'p-SH', 'c-ELW', 'p-AU']
sh_index = iqtree_topline.index("p-SH")
au_index = iqtree_topline.index("p-AU")
model_order    = ["m1", "m2", "m3", "m4", "m5", "poisson", "GTR20"]

type = sys.argv[1]
inpath = "../topology_tests_" + type + "/"
outfile = "../results_analysis/csv_files/topology_tests_" + type + ".csv"

if type == "simulation":
    outstring = "name,tree,rep,whichtest,m1,m2,m3,m4,m5,poisson,GTR20,true\n"
if type == "pandit":
    outstring = "name,whichtest,m1,m2,m3,m4,m5,poisson,GTR20\n"

#HIV_prum_rep5_AA.topology_tests          NP_spiralia_rep5_AA.topology_tests

files = [x for x in os.listdir(inpath) if x.endswith("topology_tests")]
for file in files:
    rawname = file.split(".topology_tests")[0]
    print(rawname)
    if type == "simulation":
        name = rawname.split("_AA")[0].replace("_",",").replace("rep","")
    else:
        name = rawname
    with open(inpath + file, "r") as f:
        lines = f.readlines()
    this_out = parse_tests(lines, name)
    outstring += parse_tests(lines, name)
with open(outfile, "w") as f:
    f.write(outstring)

    
