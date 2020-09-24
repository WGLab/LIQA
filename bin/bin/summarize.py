#!/usr/bin/python

from __future__ import print_function # load print function in python3
from collections import defaultdict
import math, sys, os

# set up auto dictionary function
def auto_dict():
    return defaultdict(auto_dict)


exonResultsFile = sys.argv[1]
geneResultsFile = sys.argv[2]
exonIncLvlsFile = sys.argv[3]
isoformEstsFile = sys.argv[4]
exonDASResults = sys.argv[5]
geneDASResults = sys.argv[6]

### load estimated exon inclusion levels
exonIncLvls = auto_dict()
with open(exonIncLvlsFile,"r") as FP:
    for line in FP:
        line = line.strip("\n")
        estinf = line.split("\t")

        gene = estinf[0]
        DS = estinf[1]
        cdt = estinf[2]
        est = estinf[3]
        estinf[4] = estinf[4].rstrip(";")
        exonLocations = estinf[4]

        if not bool(exonIncLvls[gene][DS][cdt]): exonIncLvls[gene][DS][cdt] = est + ","
        else: exonIncLvls[gene][DS][cdt] = exonIncLvls[gene][DS][cdt] + est + ","

        
### exon inclusion level differences calculation
diffExonIncLvls = auto_dict()
for gene in exonIncLvls:
    for DS in exonIncLvls[gene]:
        avgcdt1 = 0
        avgcdt2 = 0
        for cdt in range(1, 3):
            esttmp = exonIncLvls[gene][DS][str(cdt)].rstrip(",")
            incLvls = esttmp.split(",")
            for i in range(len(incLvls)):
                if cdt == 1: avgcdt1 += float(incLvls[i])
                if cdt == 2: avgcdt2 += float(incLvls[i])
            if cdt == 1: avgcdt1 = avgcdt1/len(incLvls)
            if cdt == 2: avgcdt2 = avgcdt2/len(incLvls)

        diff = avgcdt1 - avgcdt2
        diffExonIncLvls[gene][DS] = diff
        

### load estimated isoform relative abundances
isoEsts = auto_dict()
with open(isoformEstsFile, "r") as FP:
    for line in FP:
        line = line.strip("\n")
        estinf = line.split("\t")

        gene = estinf[0]
        isoform = estinf[1]
        est = estinf[len(estinf)-2]
        cdt = estinf[len(estinf)-1]

        if not bool(isoEsts[gene][isoform][cdt]): isoEsts[gene][isoform][cdt] = est + ","
        else: isoEsts[gene][isoform][cdt] = isoEsts[gene][isoform][cdt] + est + ","

### Hellinger distances calculation
diffHellingerDist = auto_dict()
isoAvsEsts = auto_dict()
for gene in isoEsts:
    for isoform in isoEsts[gene]:
        isoAvsEsts[gene][isoform]["1"] = 0
        isoAvsEsts[gene][isoform]["2"] = 0

        for cdt in range(1, 3):
            esttmp = isoEsts[gene][isoform][str(cdt)].rstrip(",")
            iEsts = esttmp.split(",")
            for i in range(len(iEsts)):
                isoAvsEsts[gene][isoform][str(cdt)] += float(iEsts[i])
            isoAvsEsts[gene][isoform][str(cdt)] = isoAvsEsts[gene][isoform][str(cdt)]/len(iEsts) # Average estimated isoform relative abundances

for gene in isoAvsEsts:
    diffHellingerDist[gene] = 0
    for isoform in isoAvsEsts[gene]:
        diffHellingerDist[gene] += (math.sqrt(isoAvsEsts[gene][isoform]["1"]) - math.sqrt(isoAvsEsts[gene][isoform]["2"]))**2

    diffHellingerDist[gene] = math.sqrt(diffHellingerDist[gene]/2)  # Hellinger distance differences


############################    
### OUTPUT FINAL RESULTS ###
############################

### output exon results
OUTEXON = open(exonDASResults, "w")
headerExonResults = "\t".join(["Gene","DASGroup", "ExonLocation", "Pvalue", "ExonIncLvlDiff", "Direction", "cdt1_numSamples","cdt1_ExonIncLvls", "cdt2_numSamples", "cdt2_ExonIncLvls"])
print(headerExonResults, file=OUTEXON)
with open(exonResultsFile, "r") as FP:
    for line in FP:
        line = line.strip("\n")
        testResults = line.split("\t")

        gene = testResults[0]
        DS = testResults[1]
        pvalue = testResults[2]
        testResults[3] = testResults[3].rstrip(";")
        exonLocations = testResults[3].split(";")
        diff = diffExonIncLvls[gene][DS]
        sign = "+" if diff > 0 else "-"
        diff = abs(diff)
        esttmp = exonIncLvls[gene][DS]["1"].split(",")
        numSpCdt1 = len(esttmp) - 1
        esttmp = exonIncLvls[gene][DS]["2"].split(",")
        numSpCdt2 = len(esttmp) - 1
        
        for i in range(len(exonLocations)):
            OutPut = "\t".join([gene, DS, exonLocations[i], pvalue, str(diff), sign, str(numSpCdt1) ,exonIncLvls[gene][DS]["1"], str(numSpCdt2), exonIncLvls[gene][DS]["2"]])
            print(OutPut, file=OUTEXON)                     


### output gene results
OUTGENE = open(geneDASResults, "w")
headerGeneResults = "\t".join(["Gene", "Pvalue", "HellingerDistance"])
print(headerGeneResults, file=OUTGENE)
with open(geneResultsFile, "r") as FP:
    for line in FP:
        line = line.strip("\n")
        testResults = line.split("\t")

        gene = testResults[0]
        pvalue = testResults[1]
        hellingerDist = diffHellingerDist[gene]

        OutPut = "\t".join([gene, pvalue, str(hellingerDist)])
        print(OutPut, file=OUTGENE)
