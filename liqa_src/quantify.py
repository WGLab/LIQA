#!/usr/bin/python

from __future__ import print_function # load print function in python3
from collections import defaultdict
import math, sys, os, re, pysam, time
from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter()
# set up auto dictionary function
def auto_dict():
    return defaultdict(auto_dict)


###############################################################################
###  ARGUMENT SETTINGS
###############################################################################

# checking whether argument is valid or not
validArgList = ["-bam", "-ref", "-out", "-mismatch", "-f_weight"]
for argIndex in range(1,len(sys.argv)):
    if sys.argv[argIndex][0] == "-" and sys.argv[argIndex] not in validArgList :
        print("Argument \'"+sys.argv[argIndex]+"\' is invalid!")
        sys.exit()
        

bamFileExists = 0
refFileExists = 0
outFileExists = 0
misFileExists = 0
weightFileExists = 0
for argIndex in range(1,len(sys.argv)):
    if sys.argv[argIndex] == "-bam":  ## load in BAM file
        argIndex += 1
        bamFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        bamTmp = sys.argv[argIndex].split("/")
        bamFile = bamFileAbsPath + "/" + bamTmp[len(bamTmp)-1]
        bamFileExists = 1
    elif sys.argv[argIndex] == "-ref":  ## load in annotation file
        argIndex += 1
        refFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        refTmp = sys.argv[argIndex].split("/")
        refGeneFile = refFileAbsPath + "/" + refTmp[len(refTmp)-1]
        refFileExists = 1
    elif sys.argv[argIndex] == "-out":  ## load in annotation file
        argIndex += 1
        outFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        outTmp = sys.argv[argIndex].split("/")
        outFile = outFileAbsPath + "/" + outTmp[len(outTmp)-1]
        outFileExists = 1
    elif sys.argv[argIndex] == "-mismatch":  ## mismatch tolerate
        argIndex += 1
        #misFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        #misTmp = sys.argv[argIndex].split("/")
        misMatch = int(sys.argv[argIndex])
        misFileExists = 1
    elif sys.argv[argIndex] == "-f_weight":  ## weight of F function
        argIndex += 1
        #misFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        #misTmp = sys.argv[argIndex].split("/")
        weightF = float(sys.argv[argIndex])
        weightFileExists = 1

                                                

if bamFileExists == 0 or refFileExists == 0 or outFileExists == 0 or misFileExists == 0 or weightFileExists == 0: ## lack enough arguments
    print("Please provide arguments:")
    print("-bam\tIndexed bam file")
    print("-ref\tGene annotation file")
    print("-out\tOutput file")
    print("-max_distance\tMax distance")
    print("-f_weight\tWeight of F function")
    sys.exit()


# load gene information
geneStructureInformation = auto_dict()
geneLineCount = auto_dict()

with open(refGeneFile, "r") as FP:
    for line in FP:
        line = line.strip("\n")
        tmpinf = line.split("\t")
        gene = tmpinf[0]
        
        if not bool(geneStructureInformation[gene]):
            geneLineCount[gene] = 0
            geneStructureInformation[gene][geneLineCount[gene]] = line
        else:
            geneLineCount[gene] += 1
            geneStructureInformation[gene][geneLineCount[gene]] = line

#if weightF <= 1:
    #weightF = 0
infthreshold = float(weightF)/10.0

#####################################
## Using pysam to read in bam file !!
#####################################
bamFilePysam = pysam.Samfile(bamFile,"rb")


## RESULTS FILE
OUT = open(outFile, 'w')


###########################################################################################################################
###  START TO ANALYZE DATA FOR EACH GENE ###
##########################################################################################################################

geneCount = 0

startTime = time.time()

#OUT.write("GeneName\tIsoformName\tNumberOfReads\tRelativeAbundance\n") ## Header of Results
OUT.write("GeneName\tIsoformName\tReadPerGene_corrected\tRelativeAbundance\tinfor_ratio\n")

for gene in geneStructureInformation:

    geneCount += 1
    tmpTime = (time.time() - startTime)/60.0
    
    
    sameReadCount = auto_dict()
    readStart = auto_dict()
    readEnd = auto_dict()
    readCigar = auto_dict()

    numofExons = geneLineCount[gene]
    tmpgeneinf = geneStructureInformation[gene][0].split("\t")
    geneChr = tmpgeneinf[1]
    geneStart = int(tmpgeneinf[3])
    geneEnd = int(tmpgeneinf[4])
    
    # deal with gene and isoform length information in refgene file
    tmpisoinf = tmpgeneinf[5].split(";")
    tmpgeneinf[5] = tmpisoinf[0]
    if len(tmpisoinf) == 2:
        tmpisolength = tmpisoinf[1].split(",")
        genelength = int(tmpisolength[0])-geneStart
        for iii in range(len(tmpisolength)-1):
            tmpisolength[iii] = (int(tmpisolength[iii]) - geneStart)
            if tmpisolength[iii] == 666:
                tmpisolength[iii] = 0
            else:
                tmpisolength[iii] = tmpisolength[iii]/2566
            #print(gene+"\t"+str(tmpisolength[iii]))
    else:
        genelength = 100
        tmpisolength = tmpisoinf[0].split(",")
        for iii in range(len(tmpisolength)-1):
            tmpisolength[iii] = 0

    ## load all reads information which were mapped to the specific gene within this loop using pysam
    for read in bamFilePysam.fetch(geneChr, geneStart, geneEnd):
        line = str(read)
        tmpinf = line.split("\t")
        tmpReadName = tmpinf[0]
        tmpReadChr = geneChr
        tmpReadStart = int(tmpinf[3]) + 1
        tmpReadCigar = ""

        ## Adjust to different Pysam Version!! ##

        if ")]" in tmpinf[5]: ## vector format
            
            tmpinf[5] = tmpinf[5].rstrip(")]")
            tmpinf[5] = tmpinf[5].lstrip("[(")
            tmpinfcigar = tmpinf[5].split("), (")
            for cc in tmpinfcigar:
                ttcc = cc.split(", ")
                if ttcc[0] == "3":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "N"
                if ttcc[0] == "2":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "D"
                if ttcc[0] == "1":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "I"
                if ttcc[0] == "0":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "M"
                if not (ttcc[0] == "3" or ttcc[0] == "2" or ttcc[0] == "1" or ttcc[0] == "0"):
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "X"
        else:      ## 100M10N100M format
            tmpReadCigar = tmpinf[5]
                                    
        if not bool(sameReadCount[tmpReadName]):
            sameReadCount[tmpReadName] = 1
        else:
            sameReadCount[tmpReadName] += 1
                                        
        readStart[tmpReadName][sameReadCount[tmpReadName]] = tmpReadStart
        readCigar[tmpReadName][sameReadCount[tmpReadName]] = tmpReadCigar.replace('=', "M")


    ## load structure information of the specific gene within this loop                    
                        
    tmpgeneinf[5] = tmpgeneinf[5].rstrip(",")
    isoformNames = tmpgeneinf[5].split(",")
    exonStarts = [None] * numofExons
    exonEnds = [None] * numofExons
    exonIndicators = auto_dict()
    
    for i in range(1,numofExons+1):
        tmpinf = geneStructureInformation[gene][i].split("\t")
        exonStarts[i-1] = int(tmpinf[3])+1
        exonEnds[i-1] = int(tmpinf[4])
        tmpinf[5] = tmpinf[5].rstrip(",")
        tmpExonIndicators = tmpinf[5].split(",")

        for j in range(len(tmpExonIndicators)):
            exonIndicators[isoformNames[j]][i-1] = int(tmpExonIndicators[j])

    lociIndicators = auto_dict()
    for i in range(len(isoformNames)):
        for j in range(len(exonStarts)):
            if exonIndicators[isoformNames[i]][j] == 1:
                for k in range(exonStarts[j], exonEnds[j]+1):
                    lociIndicators[isoformNames[i]][k] = 1

    #########################################################################################################################################
    ## START TO ANALYZE EACH READ 
    ##################################################################################################################################################

    qualifiedRead = auto_dict()
    prereadCount = 0
    readCount = 0
    readCount1iso = 0
    fragmentStart = auto_dict()
    fragmentEnd = auto_dict()
    CompatibleMatrix = auto_dict()
    tmpCompatibleMatrix = auto_dict()
    qualitycheck = auto_dict()
    
    for readName in sameReadCount:
        
        # load CIGAR information
        cigarNumberRead1 = auto_dict()
        cigarNumberRead2 = auto_dict()
        cigarMatchRead1 = auto_dict()
        cigarMatchRead2 = auto_dict()
        cigarInfCountRead1 = 0
        cigarInfCountRead2 = 0
        cigarInfCountRead1tmp = 0
        cigarInfCountRead2tmp = 0
        qualitycheck[readName] = 0
        
        tmp1 = re.split("([A-Z])",readCigar[readName][1])
        for i in range(len(tmp1)-1):
            if tmp1[i].isalpha():
                cigarMatchRead1[cigarInfCountRead1] = tmp1[i]
                cigarInfCountRead1 += 1
            else:
                cigarNumberRead1[cigarInfCountRead1] = int(tmp1[i])
                cigarInfCountRead1tmp += 1
                
        if sameReadCount[readName] == 2:
            tmp2 = re.split("([A-Z])",readCigar[readName][2])
            for i in range(len(tmp2)-1):
                if tmp2[i].isalpha():
                    cigarMatchRead2[cigarInfCountRead2] = tmp2[i]
                    cigarInfCountRead2 += 1
                else:
                    cigarNumberRead2[cigarInfCountRead2] = int(tmp2[i])
                    cigarInfCountRead2tmp += 1
                    
        # calculate read end positions
        readEnd[readName][1] = readStart[readName][1]
        for i in range(cigarInfCountRead1):
            readEnd[readName][1] += cigarNumberRead1[i]
            
        if sameReadCount[readName] == 2:
            readEnd[readName][2] = readStart[readName][2]
            for i in range(cigarInfCountRead2):
                readEnd[readName][2] += cigarNumberRead2[i]

        # calculate fragment START and END positions
        if sameReadCount[readName] == 2:
            fragmentStart[readName] = readStart[readName][2] if readStart[readName][1] >= readStart[readName][2] else readStart[readName][1]
            fragmentEnd[readName] = readEnd[readName][1] if readEnd[readName][1] >= readEnd[readName][2] else readEnd[readName][2]

        if sameReadCount[readName] == 1:
            fragmentStart[readName] = readStart[readName][1]
            fragmentEnd[readName] = readEnd[readName][1]
            readStart[readName][2] = readStart[readName][1]
            readEnd[readName][2] = readEnd[readName][1]
           
        ##################################################################################################################################    
        ## Obtain compatible matrix of isoforms with respect to reads
        #################################################################################################################################
        if (readStart[readName][1] >= geneStart and readStart[readName][1] <= geneEnd) and sameReadCount[readName]==1 or (readStart[readName][2] >= geneStart and readStart[readName][2] <= geneEnd and sameReadCount[readName]==2) :
            if cigarInfCountRead1 == cigarInfCountRead1tmp and cigarInfCountRead2 == cigarInfCountRead2tmp:
                base1 = readStart[readName][1] - 1
                exonIndicatorRead1 = [0] * numofExons
                if sameReadCount[readName] == 2:
                    base2 = readStart[readName][2] - 1
                    exonIndicatorRead2 = [0] * numofExons
                compatibleVector = [1] * len(isoformNames)

                ##############################################################################################################################################
                ### SET TUP COMPATIBLE INDICATOR VECTOR ###############
                ###############################################################################################################################################
                ## READ 1 ##
                # find exons where read 1 mapped to
                for i in range(cigarInfCountRead1):
                    
                    if cigarMatchRead1[i] == "M": ## matched CIGAR

                        for j in range(1,cigarNumberRead1[i]+1):
                            tmpbase = base1 + j
                            for k in range(len(exonStarts)):
                                if exonIndicatorRead1[k] > misMatch: continue
                                if tmpbase >= exonStarts[k] and tmpbase <= exonEnds[k]: exonIndicatorRead1[k] += 1 ## confirm that the read covers this exon
        
                        base1 += cigarNumberRead1[i] # jump to next match information

                    if cigarMatchRead1[i] == "N" or cigarMatchRead1[i] == "D": ## skipping area
                        base1 += cigarNumberRead1[i] # jump to next match information directly

                # set up indicator vector
                tmpcount1 = 0
                tmpcount11 = 0 ## these two variable are used to rule out skipping exons
                for i in range(len(exonIndicatorRead1)):
                    if exonIndicatorRead1[i] > misMatch: tmpcount1 += 1
                for i in range(len(exonIndicatorRead1)):

                    if exonIndicatorRead1[i] > misMatch:
                        tmpcount11 += 1
                        for j in range(len(isoformNames)):
                            if exonIndicators[isoformNames[j]][i] == 0: compatibleVector[j] = 0 ## rule out isoform j if reads covers skipping area of isoform j

                    if exonIndicatorRead1[i] <= misMatch: #aim to rule out isforms which includes exons which skipped by read
                        if tmpcount1 > 1 and tmpcount11 >= 1 and tmpcount11 < tmpcount1: ## confirm the exon i is skipped by read!!
                            for j in range(len(isoformNames)):
                                if exonIndicators[isoformNames[j]][i] == 1: compatibleVector[j] = 0

                    
                ## READ 2 ## SAME AS READ 1
                tmpcount2 = 0
                if sameReadCount[readName] == 2: ## ONLY WHEN THE READ IS PAIRED-END READ!!!
                    # find exons where read 2 mapped to
                    for i in range(cigarInfCountRead2):
                        
                        if cigarMatchRead2[i] == "M": ## matched CIGAR
                            
                            for j in range(1,cigarNumberRead2[i]+1):
                                tmpbase = base2 + j
                                for k in range(len(exonStarts)):
                                    if exonIndicatorRead2[k] > misMatch: continue
                                    if tmpbase >= exonStarts[k] and tmpbase <= exonEnds[k]: exonIndicatorRead2[k] += 1 ## confirm that the read covers this exon
                                    
                            base2 += cigarNumberRead2[i] # jump to next match information
                                    
                        if cigarMatchRead2[i] == "N" or cigarMatchRead1[i] == "D": ## skipping area
                            base2 += cigarNumberRead2[i] # jump to next match information directly
                                        
                    # set up indicator vector
                    tmpcount2 = 0
                    tmpcount22 = 0 ## these two variable are used to rule out skipping exons
                    for i in range(len(exonIndicatorRead2)):
                        if exonIndicatorRead2[i] > misMatch: tmpcount2 += 1
                    for i in range(len(exonIndicatorRead2)):
                            
                        if exonIndicatorRead2[i] > misMatch:
                            tmpcount22 += 1
                            for j in range(len(isoformNames)):
                                if exonIndicators[isoformNames[j]][i] == 0: compatibleVector[j] = 0 ## rule out isoform j if reads covers skipping area of isoform j
                                                    
                        if exonIndicatorRead2[i] <= misMatch: #aim to rule out isforms which includes exons which skipped by read
                            if tmpcount2 > 1 and tmpcount22 >= 1 and tmpcount22 < tmpcount2: ## confirm the exon i is skipped by read!!
                                for j in range(len(isoformNames)):
                                    if exonIndicators[isoformNames[j]][i] == 1: compatibleVector[j] = 0
                                                                
                ##################################################################################################################################################
                ## fill in compatible matrix ##
                if tmpcount1 > 0 or (tmpcount2 > 0 and sameReadCount[readName] == 2):
                    prereadCount += 1
                    for i in range(len(isoformNames)):
                        CompatibleMatrix[readName][isoformNames[i]] = compatibleVector[i]
                        tmpCompatibleMatrix[readName][isoformNames[i]] = compatibleVector[i]
                        qualitycheck[readName] += compatibleVector[i]
                    if qualitycheck[readName] > 0:
                        qualifiedRead[readName] = 1
                        #readCount += 1
                    if qualitycheck[readName] == 1:
                        readCount1iso += 1
                    readCount += 1

    ### COMPATIBLE MATRIX OBTAINED !!!
    ###############################################################################################################
    
    if readCount == 0: continue
    fisherinf = readCount1iso/float(prereadCount)
    if fisherinf < infthreshold: continue
    print(gene+"\tprocessing...")
    
    ##############################################################################################################
    ### ANALYZE EMPIRICAL READS DISTRIBUTION BASED ON NON-PARAMATRIC METHOD
    ##############################################################################################################

    positionsFragmentCovered = auto_dict()
    readsDistributionIsoform = auto_dict()
    readsDistributionIsoformKnown = auto_dict()
    denominatorKnown = auto_dict()
    denominator = auto_dict()
    
    for i in range(len(isoformNames)):
        for readName in qualifiedRead:
            ####################################################################################################
            ### CACULATE VALID FRAGMENT LOCATION ON THIS GENE ###################################
            tmpStart = None
            tmpEnd = None
            #print(readName)
            for j in range(len(exonStarts)):
                if j == 0:
                    if fragmentStart[readName] < exonStarts[j]:
                        tmpStart = exonStarts[j]
                    if fragmentStart[readName] >= exonStarts[j] and fragmentStart[readName] <= exonEnds[j]:
                        tmpStart = fragmentStart[readName]
                    if fragmentEnd[readName] >= exonStarts[j] and fragmentEnd[readName] <= exonEnds[j]:
                        tmpEnd = fragmentEnd[readName]

                if j > 0 and j < len(exonStarts)-1:
                    if fragmentStart[readName] < exonStarts[j] and fragmentStart[readName] > exonEnds[j-1]:
                        tmpStart = exonStarts[j]
                    if fragmentStart[readName] >= exonStarts[j] and fragmentStart[readName] <= exonEnds[j]:
                        tmpStart = fragmentStart[readName]
                    if fragmentEnd[readName] < exonStarts[j] and fragmentEnd[readName] > exonEnds[j-1]:
                        tmpEnd = exonEnds[j-1]
                    if fragmentEnd[readName] >= exonStarts[j] and fragmentEnd[readName] <= exonEnds[j]:
                        tmpEnd = fragmentEnd[readName]

                if j == len(exonStarts)-1:
                    if fragmentStart[readName] >= exonStarts[j] and fragmentStart[readName] <= exonEnds[j]:
                        tmpStart = fragmentStart[readName]
                    if fragmentStart[readName] < exonStarts[j] and fragmentStart[readName] > exonEnds[j-1]:
                        tmpStart = exonStarts[j]
                    if fragmentEnd[readName] >= exonStarts[j] and fragmentEnd[readName] <= exonEnds[j]:
                        tmpEnd = fragmentEnd[readName]
                    if fragmentEnd[readName] < exonStarts[j] and fragmentEnd[readName] > exonEnds[j-1]:
                        tmpEnd = exonEnds[j-1]
                    if fragmentEnd[readName] > exonEnds[j]:
                        tmpEnd = exonEnds[j]

            ## Valid position obtained
            #################################################################

            fragmentStart[readName] = tmpStart ## new starts and end position updated
            fragmentEnd[readName] = tmpEnd

            
    #####################################################################################################
    ## EM algorithm
    #####################################################################################################

    Alpha = [1.0] * len(isoformNames)
    oldAlpha = [None] * len(isoformNames)

    isoformLength = auto_dict()
    
    for i in range(len(isoformNames)):
        for j in range(len(exonStarts)):
            if exonIndicators[isoformNames[i]][j] == 1:
                if not bool(isoformLength[isoformNames[i]]): isoformLength[isoformNames[i]] = exonEnds[j] - exonStarts[j] + 1
                else: isoformLength[isoformNames[i]] += exonEnds[j] - exonStarts[j] + 1  ## calculate isoform length

    ########################################################################################################            
    ## UPDATE empirical probability of read mapped to the specific position
 
    Hfunction = auto_dict()
    numerator = auto_dict()
    numeratorKnown = auto_dict()
    for i in range(len(isoformNames)):
        if not bool(isoformLength[isoformNames[i]]): isoformLength[isoformNames[i]] = 0
        for readName in qualifiedRead:
            if isoformLength[isoformNames[i]] > 0: Hfunction[readName][isoformNames[i]] = float(1)
            else: Hfunction[readName][isoformNames[i]] = 0

    #########################################################################################################
    ## iteration begins        
            
    diff = 1.0
    iterCount = 0
    while diff > .0001:
        #print(gene+"\t"+str(geneCount)+"\t"+str(iterCount)+"\t"+str(diff)+"\t"+str(tmpTime))
        up = None
        down = None
        if iterCount == 0:

            for readName in qualifiedRead:
                    
                    #up = 0.0
                down = 0.0
                for j in range(len(isoformNames)):
                    if CompatibleMatrix[readName][isoformNames[j]] == 1:
                        if (isoformLength[isoformNames[j]]-300+1) > 0: down += Alpha[j] / (isoformLength[isoformNames[j]]-300+1)
                        else: down += Alpha[j]
                for i in range(len(isoformNames)):
                    up = 0.0
                    if CompatibleMatrix[readName][isoformNames[i]] == 1:
                        if (isoformLength[isoformNames[i]]-300+1) > 0: up = Alpha[i] / (isoformLength[isoformNames[i]]-300+1)
                        else: up = Alpha[i]
                    if down != 0: tmpCompatibleMatrix[readName][isoformNames[i]] = up / down

        if iterCount > 0:
            for readName in qualifiedRead:

                    #up = 0.0
                down = 0.0
                for j in range(len(isoformNames)):
                    if CompatibleMatrix[readName][isoformNames[j]] == 1: down += Alpha[j] * Hfunction[readName][isoformNames[j]]
                for i in range(len(isoformNames)):
                    up = 0.0
                    if CompatibleMatrix[readName][isoformNames[i]] == 1: up = Alpha[i] * Hfunction[readName][isoformNames[i]]
                    if down != 0: tmpCompatibleMatrix[readName][isoformNames[i]] = up / down

        for i in range(len(isoformNames)):
            up1 = 0.0
            for readName in qualifiedRead:
                up1 += tmpCompatibleMatrix[readName][isoformNames[i]]
            oldAlpha[i] = Alpha[i]
            Alpha[i] = up1 / readCount

        iterCount += 1
        for t in range(len(Alpha)):
            if t == 0: diff = abs(Alpha[t] - oldAlpha[t])
            else:
                if abs(Alpha[t] - oldAlpha[t]) > diff: diff = abs(Alpha[t] - oldAlpha[t])

                
    sumAlpha = sum(Alpha)
    if sumAlpha == 0: continue

    for i in range(len(Alpha)):
        Alpha[i] = Alpha[i] / sumAlpha

    isoformRelativeAbundances = [None] * len(isoformNames)
    sumTheta = 0.0
    for i in range(len(Alpha)):
        #sumTheta += Alpha[i] / isoformLength[isoformNames[i]] + weightF * tmpisolength[i]
        sumTheta += Alpha[i] #+ weightF * tmpisolength[i]

    print(gene+"\t"+str(iterCount)+" iterations\tDone!")
    
    #rpg_lengthcorrected = readCount/genelength*100
    for i in range(len(Alpha)):
        #isoformRelativeAbundances[i] = (Alpha[i]/isoformLength[isoformNames[i]]+tmpisolength[i]*weightF) /(sumTheta)
        #isoformRelativeAbundances[i] = (Alpha[i] + tmpisolength[i]*weightF) /(sumTheta)
        isoformRelativeAbundances[i] = (Alpha[i]) /(sumTheta)
        #print(gene+"\t"+str(geneCount)+"\t"+str(iterCount)+"\t"+isoformNames[i]+"\t"+str(isoformRelativeAbundances[i])+"\t"+str(tmpTime))
        
        rpg_lengthcorrected = readCount/genelength*100*isoformRelativeAbundances[i]

        #OUT.write(gene+"\t"+isoformNames[i]+"\t"+str(readCount)+"\t"+str(isoformRelativeAbundances[i])+"\t"+str(rpg_lengthcorrected)+"\n") ## write results into specified file
        OUT.write(gene+"\t"+isoformNames[i]+"\t"+str(rpg_lengthcorrected)+"\t"+str(isoformRelativeAbundances[i])+"\t"+str(fisherinf)+"\n")

OUT.close()
            
                            


                    
