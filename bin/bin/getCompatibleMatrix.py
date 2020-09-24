#!/usr/bin/python

from __future__ import print_function # load print function in python3
from collections import defaultdict  
import sys, os

# set up auto dictionary function
def auto_dict():
    return defaultdict(auto_dict)


refGeneFile = sys.argv[1]

########################################
##### load transcript annotations ######
########################################
knownGene = auto_dict()
with open(refGeneFile,"r") as FP:
    for line in FP:

        line = line.strip("\n")
        transcript = line.split("\t")
        if "NR" in transcript[1]:
            continue       

        gene = transcript[12]
        tran = transcript[1]
        knownGene[gene][tran]["chrom"] = transcript[2]
        knownGene[gene][tran]["strand"] = transcript[3]
        knownGene[gene][tran]["txStart"] = transcript[4]
        knownGene[gene][tran]["txEnd"] = transcript[5]
        knownGene[gene][tran]["exonCount"] = transcript[8]
        knownGene[gene][tran]["exonStarts"] = transcript[9].rstrip(",")
        knownGene[gene][tran]["exonEnds"] = transcript[10].rstrip(",")

####################################
##### load isoform annotations #####
####################################
isoGene = auto_dict()
geneStart = auto_dict()
geneEnd = auto_dict()


for gene in knownGene:
    for tran in knownGene[gene]:

        g_start = knownGene[gene][tran]["txStart"]
        g_end = knownGene[gene][tran]["txEnd"]

        if not gene in isoGene:
            isoGene[gene] = tran + ","
        else:
            isoGene[gene] = isoGene[gene] + tran + ","

        if not gene in geneStart:
            geneStart[gene] = g_start
        else:
            if int(geneStart[gene]) > int(g_start): geneStart[gene] = g_start

        if not gene in geneEnd:
            geneEnd[gene] = g_end
        else:
            if int(geneEnd[gene]) < int(g_end): geneEnd[gene] = g_end


#######################################
##### process isoform information #####
#######################################

for gene in isoGene:

    isoName = isoGene[gene].split(",")
    size = len(isoName)-1
    g_chrom = knownGene[gene][isoName[0]]["chrom"]
    g_strand = knownGene[gene][isoName[0]]["strand"]
    g_start = geneStart[gene]
    g_end = geneEnd[gene]
    ISO_INDEX = auto_dict()

    if not size >= 2: ## exclude genes with # isoforms less than 2
        continue
    

    tmp = (gene,g_chrom,g_strand,g_start,g_end)
    print("\t".join(tmp)+"\t",end="")

    for j in range(0, len(isoName)-1):

        tran = isoName[j]
        print(tran+",",end="")

        isoStarts = knownGene[gene][tran]["exonStarts"].split(",")
        isoEnds = knownGene[gene][tran]["exonEnds"].split(",")

        for ijk in range(0, len(isoStarts)):
            sss = isoStarts[ijk]
            eee = isoEnds[ijk]
            for abc in range(int(sss),int(eee)+1):
                ISO_INDEX[abc][j] = 1
    
    print("\n",end="")

    ####################################
    ### get virtual exon information ###
    ####################################
    NEW_EXON = auto_dict()
    ccc = 0
    pre_pos = int(g_start) - 10
    pre_index = auto_dict()
    for j in range(0, len(isoName)-1):
        pre_index[j] = 0

    ISO_INDEX_keys_sorted = ISO_INDEX.keys()
    ISO_INDEX_keys_sorted.sort()
    for ijk in ISO_INDEX_keys_sorted:

        tot = 0
        cur_index = auto_dict()
        for j in range(0, len(isoName)-1):

            iso = isoName[j]
            if bool(ISO_INDEX[ijk][j]):
                cur_index[j] = 1
            if not bool(ISO_INDEX[ijk][j]):
                cur_index[j] = 0
            if cur_index[j] != pre_index[j]:
                tot = tot + 1
        move = ijk - pre_pos

        if move != 1: ## jump to new exon across intron
            NEW_EXON[ccc]["start"] = ijk
            NEW_EXON[ccc]["index"] = cur_index
            if ccc > 0:
                NEW_EXON[ccc-1]["end"] = pre_pos
            pre_index = cur_index
            ccc = ccc + 1
        else:
            if tot > 0: ## new virtual exon
                NEW_EXON[ccc]["start"] = ijk
                NEW_EXON[ccc-1]["end"] = ijk-1
                NEW_EXON[ccc]["index"] = cur_index
                pre_index = cur_index
                ccc = ccc + 1

        pre_pos = ijk

    ##ijk    
    NEW_EXON[ccc-1]["end"] = g_end

    ######################
    ### print out data ###
    ######################
    NEW_EXON_keys_sorted = NEW_EXON.keys()
    NEW_EXON_keys_sorted.sort()
    for ccc in NEW_EXON_keys_sorted:

        tmp = (gene,g_chrom,g_strand)
        print("\t".join(tmp)+"\t",end="")
        sss = NEW_EXON[ccc]["start"]
        print(str(sss)+"\t",end="")
        eee = NEW_EXON[ccc]["end"]
        print(str(eee)+"\t",end="")

        index_keys_sorted = NEW_EXON[ccc]["index"].keys()
        index_keys_sorted.sort()
        for k in index_keys_sorted:
            index_get = str(NEW_EXON[ccc]["index"][k])
            print(index_get+",",end="")
        print("\n",end="")














        
