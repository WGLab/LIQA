#!/usr/bin/python

from liqa_src import my_functions as my

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
        validArgList = ["-task", "-ref", "-format", "-out"]
        addAbsPath = [0, 1, 0, 3]
        message = "liqa -task refgene -ref <reference_file> -format <reference_file_format(gtf/ucsc)> -out <output_file>"
        inputs = my.parse_argument(validArgList, addAbsPath, message)
        refFile = inputs[1]
        formatFile = inputs[2]
        outFile = inputs[3]

        if (formatFile == "ucsc"):
            myCommand = "perl " + fileAbsPath + "/PreProcess.pl -r " + refFile + " -o " + outFile
            os.system(myCommand)
        elif (formatFile == "gtf"):
            myCommand = "perl " + fileAbsPath + "/PreProcess_gtf.pl -r " + refFile + " -o " + outFile
            os.system(myCommand)
        else:
            print("Please specify reference file format: gtf/ucsc")

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
        myCommand = "python " + fileAbsPath + "/quantify.py -ref " + refFile + " -bam " +  bamFile + " -out " + outFile + " -mismatch " + misMatch + " -f_weight " + weightF
        #myCommand = "python " + "liqa_bin/quantify.py -ref " + refFile + " -bam " +  bamFile + " -out " + outFile + " -mismatch " + misMatch + " -f_weight " + weightF
        os.system(myCommand)

    if task == "diff":
        validArgList = ["-task", "-condition_1", "-condition_2", "-out"]
        addAbsPath = [0, 1, 1, 3]
        message = "liqa    -task diff\n\t-condition_1 <isoform_expression_estimation_file_for_condition1>\n\t-condition_2 <isoform_expression_estimation_file_for_condition2>\n\t-out <test_results_file>"
        inputs = my.parse_argument(validArgList, addAbsPath, message)
        cdt1 = inputs[1]
        cdt2 = inputs[2]
        outFile = inputs[3]
        myCommand = "perl " + fileAbsPath + "/group_process.pl -gp1 " + cdt1 + " -gp2 " + cdt2 + " -o " + crtAbsPath + "/isoform_expression_summary"
        #myCommand = "perl " + "liqa_bin/group_process.pl -gp1 " + cdt1 + " -gp2 " + cdt2 + " -o " + crtAbsPath + "/isoform_expression_summary"
        os.system(myCommand)
        
        myCommand = "Rscript " + fileAbsPath + "/testDAS.R " + crtAbsPath + "/isoform_expression_summary " + outFile
        #myCommand = "Rscript " + "liqa_bin/testDAS.R " + crtAbsPath + "/isoform_expression_summary " + outFile
        os.system(myCommand)

if __name__ == "__main__":
    main()

