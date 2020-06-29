#!/usr/bin/python

from bin import my_functions as my

import os, sys
fileAbsPath = os.path.abspath(os.path.dirname(__file__))
crtAbsPath = os.getcwd()

task = ""
taskList = ["refgene", "quantify"]
for i in range(1,len(sys.argv)):
    if sys.argv[i] == "-task" and len(sys.argv)!=i+1:
        task = sys.argv[i+1]
if (task not in taskList):
    print("\nPlease specify task (LIQA.py -task <task>):\n")
    print("\trefgene:   preprocess reference file\n")
    print("\tquantify:   quantify isoform expression\n")

if task == "refgene":
    validArgList = ["-task", "-ref", "-out"]
    addAbsPath = [0, 1, 3]
    message = "LIQA.py -task refgene -ref <reference_file> -out <output_file>"
    inputs = my.parse_argument(validArgList, addAbsPath, message)
    refFile = inputs[1]
    outFile = inputs[2]
    myCommand = "perl " + fileAbsPath + "/bin/PreProcess.pl -r " + refFile + " -o " + outFile
    os.system(myCommand)

if task == "quantify":
    validArgList = ["-task", "-refgene", "-bam", "-out"]
    addAbsPath = [0, 1, 1, 3]
    message = "LIQA.py -task quantify -refgene <refgene_file> -bam <bam_file> -out <output_file>"
    inputs = my.parse_argument(validArgList, addAbsPath, message)
    refFile = inputs[1]
    bamFile = inputs[2]
    outFile = inputs[3]
    myCommand = "python " + fileAbsPath + "/bin/LRSeq_new.py -ref " + refFile + " -bam " +  bamFile + " -out " + outFile
    os.system(myCommand)
