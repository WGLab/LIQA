#!/usr/bin/python

from liqa_bin import my_functions as my

import os, sys

def main():
    
    fileAbsPath = os.path.abspath(os.path.dirname(__file__))
    crtAbsPath = os.getcwd()

    task = ""
    taskList = ["refgene", "quantify", "diff"]
    for i in range(1,len(sys.argv)):
        if sys.argv[i] == "-task" and len(sys.argv)!=i+1:
            task = sys.argv[i+1]
    if (task not in taskList):
        print("\nPlease specify task (liqa -task <task>):\n")
        print("\trefgene:   preprocess reference file\n")
        print("\tquantify:   quantify isoform expression\n")
        print("\tdiff:   detect differential splicing gene/isoform\n")

    if task == "refgene":
        validArgList = ["-task", "-ref", "-out"]
        addAbsPath = [0, 1, 3]
        message = "liqa -task refgene -ref <reference_file> -out <output_file>"
        inputs = my.parse_argument(validArgList, addAbsPath, message)
        refFile = inputs[1]
        outFile = inputs[2]
        myCommand = "perl " + fileAbsPath + "/liqa_function/PreProcess.pl -r " + refFile + " -o " + outFile
        os.system(myCommand)

    if task == "quantify":
        validArgList = ["-task", "-refgene", "-bam", "-out", "-max_distance", "-f_weight"]
        addAbsPath = [0, 1, 1, 3, 0, 0]
        message = "liqa -task quantify -refgene <refgene_file> -bam <bam_file> -out <output_file> -max_distance <max distance> -f_weight <weight of F function>"
        inputs = my.parse_argument(validArgList, addAbsPath, message)
        refFile = inputs[1]
        bamFile = inputs[2]
        outFile = inputs[3]
        misMatch = inputs[4]
        weightF = inputs[5]
        myCommand = "python " + fileAbsPath + "/liqa_function/LRSeq_new.py -ref " + refFile + " -bam " +  bamFile + " -out " + outFile + " -mismatch " + misMatch + " -f_weight " + weightF
        os.system(myCommand)

    if task == "diff":
        validArgList = ["-task", "-ref", "-est"]
        addAbsPath = [0, 1, 1]
        message = "liqa -task refgene -ref <refgene_file> -est <isoformRelativeAbundances_estimations>"
        inputs = my.parse_argument(validArgList, addAbsPath, message)
        refFile = inputs[1]
        estFile = inputs[2]
        myCommand = "python " + fileAbsPath + "/liqa_function/Diff.py -r " + refFile + " -est " + estFile
        os.system(myCommand)

if __name__ == "__main__":
    main()

