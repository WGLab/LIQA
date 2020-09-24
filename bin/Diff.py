#!/usr/bin/python

from __future__ import print_function # load print function in python3
from collections import defaultdict  
import sys, os

# set up auto dictionary function
def auto_dict():
    return defaultdict(auto_dict)

# preliminary setup
scriptAbsPath = os.path.dirname(os.path.abspath(__file__)) # get absolute path of python script
while not os.path.isdir(scriptAbsPath+"/tmp"): os.system("mkdir "+scriptAbsPath+"/tmp")

# checking out argments
validArgList = ["-est", "-ref"]
for argIndex in range(1,len(sys.argv)):
    if sys.argv[argIndex][0] == "-" and sys.argv[argIndex] not in validArgList :
        print("Argument \'"+sys.argv[argIndex]+"\' is invalid!")
        sys.exit()
        
estFileExists = 0
refFileExists = 0
for argIndex in range(1,len(sys.argv)):
    if sys.argv[argIndex] == "-est":  ## load in estimation file
        argIndex += 1
        estFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        estTmp = sys.argv[argIndex].split("/")
        estimationFile = estFileAbsPath + "/" + estTmp[len(estTmp)-1]
        estFileExists = 1
    elif sys.argv[argIndex] == "-ref":  ## load in annotation file
        argIndex += 1
        refFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        refTmp = sys.argv[argIndex].split("/")
        refGeneFile = refFileAbsPath + "/" + refTmp[len(refTmp)-1]
        refFileExists = 1

if estFileExists == 0 or refFileExists == 0: ## lack enough arguments
    print("Please provide arguments:")
    print("-est\tIsoform relative abundances estimation file")
    print("-ref\tGene annotation file")
    sys.exit()    



### convert gene annotation file to compatible matrix

getRefGeneCommand = "python " + scriptAbsPath + "/bin/getCompatibleMatrix.py " + refGeneFile + " > " + scriptAbsPath + "/tmp/isoformCompatibleMatrix"  
command1Done = 0
while command1Done == 0:
    print("Isoform compatible matrix conversion: Processing...")
    os.system(getRefGeneCommand) 
    command1Done = 1

print("Isoform compatible matrix conversion: Done!")


### calculate estimated exon inclusion levels
getExonIncLvlCommand = "python " + scriptAbsPath + "/bin/getExonIncLvl.py " + scriptAbsPath + "/tmp/isoformCompatibleMatrix" + " " + estimationFile + " > " + scriptAbsPath + "/tmp/exonInclusionLevels"
command2Done = 0
while command2Done == 0:
    print("Exon inclusion levels calculation: Processing...")
    os.system(getExonIncLvlCommand)
    command2Done = 1

print("Exon inclusion levels calculation: Done!")

### test DAS event
testDASCommand = "Rscript " + scriptAbsPath + "/bin/testDAS.R " + scriptAbsPath + "/tmp/exonInclusionLevels " + scriptAbsPath + "/tmp/tmpResultsGene "+ scriptAbsPath +"/tmp/tmpResultsExon"
command3Done = 0
while command3Done == 0:
    print("DAS testing: Start to process...")
    os.system(testDASCommand)
    command3Done = 1

print("DAS testing: Done!")


### summarize results
summarizeCommand = "python " + scriptAbsPath + "/bin/summarize.py " + scriptAbsPath + "/tmp/tmpResultsExon " + scriptAbsPath + "/tmp/tmpResultsGene "
summarizeCommand = summarizeCommand + scriptAbsPath + "/tmp/exonInclusionLevels " + estimationFile + " " + scriptAbsPath + "/exonDASResults " + scriptAbsPath + "/geneDASResults"
command4Done = 0
while command4Done == 0:
    print("Summarizing...")
    os.system(summarizeCommand)
    command4Done = 1

print("Done!")

### delete temporary files
deleteCommand = "rm -r " + scriptAbsPath + "/tmp"
os.system(deleteCommand)
