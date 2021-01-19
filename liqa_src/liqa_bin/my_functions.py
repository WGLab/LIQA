#!/usr/bin/python

from __future__ import print_function
from collections import defaultdict
import math, sys, os, re, time

# set up auto dictionary function
def auto_dict():
    return defaultdict(auto_dict)

# make a directory
def mk_dir(path):
    check = os.path.isdir(path)
    if not check:
        os.system("mkdir " + path)
    return

# parse arguments
def parse_argument(validArgList, addAbsPath, warnMessage):    
    for argIndex in range(1,len(sys.argv)):
        if sys.argv[argIndex][0] == "-" and sys.argv[argIndex] not in validArgList :
            print("Argument \'"+sys.argv[argIndex]+"\' is invalid!")
            sys.exit()

    # assign arguments to a list
    outList = []
    for i in range(0, len(validArgList)):
        for argIndex in range(1,len(sys.argv)):
            if sys.argv[argIndex] == validArgList[i]:
                argIndex += 1
                if "~" in sys.argv[argIndex]:
                    sys.argv[argIndex] = os.path.expanduser(sys.argv[argIndex])
                fileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
                fileTmp = sys.argv[argIndex].split("/")
                if addAbsPath[i] == 1: # target file
                    fileTmp = fileAbsPath + "/" + fileTmp[len(fileTmp)-1]
                    check = os.path.exists(fileTmp)
                    if not check:
                        print(fileTmp+" does not exist!")
                        sys.exit()
                if addAbsPath[i] == 3: # create target file
                    fileTmp = fileAbsPath + "/" + fileTmp[len(fileTmp)-1]
                if addAbsPath[i] == 0: # value
                    fileTmp = fileTmp[len(fileTmp)-1]
                if addAbsPath[i] == 2: # target directory
                    fileTmp = os.path.abspath(sys.argv[argIndex])
                    check = os.path.isdir(fileTmp)
                    if not check:
                        print(fileTmp+" does not exist!")
                        sys.exit()
                outList.append(fileTmp)
    
    if len(outList) != len(validArgList):
        print(warnMessage)
        sys.exit()
    return outList

# check modules ### NOT WORKING!!
import imp
def check_module_exists(name):
    try:
        imp.find_module(name)
    except ImportError:
        return False
    return True

def check_module(module):
    x = check_module_exists(module)
    if x:
        print("Module \'" + module + "\' is installed.")
    if not x:
        print("Module \'" + module + "\' is NOT installed!")
    return

# check program
from subprocess import Popen, PIPE

def check_program_exists(name):
    p = Popen(['/usr/bin/which', name], stdout=PIPE, stderr=PIPE)
    p.communicate()
    return p.returncode == 0

def check_program(program):
    x = check_program_exists(program)
    if x:
        print("Program \'" + program + "\' is installed.")
    if not x:
        print("Program \'" + program + "\' is NOT installed!")
    return
