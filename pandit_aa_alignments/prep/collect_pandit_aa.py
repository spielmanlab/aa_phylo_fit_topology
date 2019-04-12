"""
SJS
Script to browse pandit and retrieve alignments used to derive simulation parameters.

    Pandit alignments are only considered if they meet criteria of: 
       a) >=100 sequences
       b) <=0.7 average pairwise identity
       c) No stop codons
       d) columns have >=80% data (<20% missing culled columns)
    If the final alignment contains over 750 columns (aka 250 codon positions), retain.
    If the final alignment contains more than 500 sequences, randomly downsample for swmutsel inference runtime ease.
"""

from Bio import AlignIO, SeqIO
import numpy as np
import subprocess
from urllib.request import urlopen

import re
import os

base_search_url = "http://www.ebi.ac.uk/goldman-srv/pandit/pandit.cgi?action=browse&key=REPLACE"
base_aln_url     = "https://www.ebi.ac.uk/goldman-srv/pandit/pandit.cgi?action=browse&fam=REPLACE&field=nam-dsq-asq"
outdir = "pandit_aa_alignments/"
os.system("mkdir -p "+outdir)


def obtain_source(search_url, type):

    response = urlopen(search_url)
    if type == "lines":
        page_source = [x.decode(response.headers.get_content_charset()) for x in response.readlines()]
    elif type == "full":
        page_source = response.read().decode(response.headers.get_content_charset())
        
    return page_source



def parse_pandit_page(query_pf):

    id_list = []
    
    # download the page source
    url = base_search_url.replace("REPLACE",query_pf)
    page_source = obtain_source(url, "lines")
    
    # parse the page source
    for line in page_source:
        if line.startswith('<tr bgcolor="#eeeeee">'):
            split_line = line.split('</td><td align="right">')[1:]
            numseq = int(split_line[0])
            pairid = float(split_line[2])
            find_pfam = re.search("family\?entry=(PF\d+)", split_line[-1])
            if find_pfam:
                pfam_id = find_pfam.group(1)
                id_list.append(pfam_id)
    return id_list


def main():

    for i in range(76,82):
        print( i)
        # Generate search key
        query_pf = str(i)
        while len(query_pf) != 3:
            query_pf = "0" + query_pf
        query_pf = "PF" + query_pf

        # download and parse the page source to obtain IDs we want
        ids = parse_pandit_page(query_pf)
        
        # grab the alignment and tree for each entry in ids
        for id in ids:
            if not os.path.exists(outdir + id + ".fasta"):
                url = base_aln_url.replace("REPLACE",id)
                page_source = obtain_source(url, "full")
        
                # Grab just the table lines
                page_source = page_source.split('<table cellspacing="0" cellpadding="4" width="100%" border="0">')[1].split('</table>')[0]
        
                # split table into key, value 
                page_source = page_source.replace("<tr>\n", "").replace("</tt></td>\n</tr>", "").replace('<td align="left" bgcolor="#eeeeee" nowrap><tt>', "").replace("</tt></td>", ";")
                page_source = page_source.replace("/", "_").replace(";\n", ";")
                pairs = page_source.split("\n")[1:-1]
        
                # Save to TEMP fasta file
                records = {}
                goodid = ""
                for entry in pairs:
                    idseq = entry.split(";")
                    thisid = idseq[0].replace("-","_")
                    seq = idseq[1].replace(".", "-").upper()
                    records[thisid] = seq
                    goodid = thisid

                with open(outdir + id + ".fasta", "w") as f:
                    for r in records:
                        f.write(">" + r + "\n" + records[r] + "\n")
main()
    
    

    

    
