#!/usr/bin/python

from __future__ import print_function # load print function in python3
from collections import defaultdict  
import os, sys

# set up auto dictionary function
def auto_dict():
    return defaultdict(auto_dict)



refGeneFile = sys.argv[1]
estimationFile = sys.argv[2]


########################################
##### load transcript annotations ######
########################################
DASCount = auto_dict()
DASCombination = auto_dict()
DASExons = auto_dict()
geneToTrans = auto_dict()
tmpcount = auto_dict()

with open(refGeneFile,"r") as FP:
    for line in FP:
        line = line.strip("\n")
        refInf = line.split("\t")
        
        gene = refInf[0]
        chrom = refInf[1]
        start = refInf[3]
        end = refInf[4]
        isoInf = refInf[len(refInf)-1]
        
        isoInf = isoInf.rstrip(",")
        tmpcount[gene] = 1 if not bool(tmpcount[gene]) else tmpcount[gene] + 1

        ##  load isoform names
        if tmpcount[gene] == 1:
            trans = isoInf.split(",")
            for i in range(len(trans)):
                geneToTrans[gene][i] = trans[i]                
            continue

        ## load isoform indexes
        indexes = isoInf.split(",")
        weight = 0
        if "0," in isoInf:
            exonLocation = (",").join((chrom, start, end))+";"
            tmp = ""
            tmpweight = 0
            for i in range(len(indexes)):         
                weight = weight + 2**i             ### assign weight to exon
                if indexes[i] == "1":
                    tmpweight = tmpweight + 2**i
                    tmp = tmp + str(i) + ","

            if tmpweight > weight/2.0:
                tmpweight = weight - tmpweight

            DASExons[gene][tmpweight] = exonLocation if not bool(DASExons[gene][tmpweight]) else DASExons[gene][tmpweight] + exonLocation
            DASCount[gene][tmpweight] = 1 if not bool(DASCount[gene][tmpweight]) else DASCount[gene][tmpweight] + 1    
            DASCombination[gene][tmpweight] = tmp


            
###################################
##### load estimation results #####
###################################            
estimationResults = auto_dict()
geneToConditions = auto_dict()
with open(estimationFile, "r") as FP:
    for line in FP:
        line = line.strip("\n")
        results = line.split("\t")        
        estimationResults[results[0]][results[len(results)-3]][results[1]] = results[len(results)-2]
        if results[1] == geneToTrans[results[0]][0]:
            geneToConditions[gene] = results[len(results)-1] + "," if not bool(geneToConditions[gene]) else geneToConditions[gene]+results[len(results)-1] + ","


######################################################
##### calculate exon inclusion level estimations #####
######################################################
indexCount = auto_dict()
for gene in estimationResults:
    eventnum = 0
    geneToConditions[gene] = geneToConditions[gene].rstrip(",")
    conditions = geneToConditions[gene].split(",")

    for weight in DASCount[gene]:
        eventnum = eventnum + 1
        DASCombination[gene][weight] = DASCombination[gene][weight].rstrip(",")
        trans = DASCombination[gene][weight].split(",")

        for sp in range(len(conditions)):
            if bool(estimationResults[gene][str(sp)]):
                tmpest = 0

                for i in range(len(trans)):
                    tmpest = tmpest + float( estimationResults[gene][str(sp)][geneToTrans[gene][int(trans[i])]]  )


                indexCount[gene][eventnum] = 1 if not bool(indexCount[gene][eventnum]) else indexCount[gene][eventnum] + 1
                outPut = (gene, "DS"+str(eventnum), conditions[sp], str(tmpest), DASExons[gene][weight], str(indexCount[gene][eventnum])  )
                print("\t".join(outPut))
