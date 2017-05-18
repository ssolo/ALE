#!/usr/bin/python
import sys
import subprocess

#print(sys.argv)
#print (sys.argv[1:len(sys.argv)])

prog_names = list()
prog_names.append("ALEobserve")
prog_names.append("ALEml")
#prog_names.append("ALEmcmc")
prog_names.append("ALEml_undated")
prog_names.append("ALEmcmc_undated")
prog_names.append("ALEsample")
prog_names.append("ALEadd")
prog_names.append("ALEcount")
prog_names.append("CCPscore")
prog_names.append("computeALEcomplexity")
prog_names.append("ls_leaves")
prog_names.append("mlresampler_undated")
prog_names.append("simulation")
prog_names.append("times")


if (len(sys.argv)==1):
    print("Problem, you need to provide an ALE program name.\n")
    print("Admissible program names are: "+str(prog_names) +"\n\n")
elif (sys.argv[1] in prog_names):
    print("\n\t\tRunning "+ sys.argv[1] + ":\n")
    sys.argv[1] = "/usr/local/ALE/build/bin/"+sys.argv[1]
    subprocess.call(sys.argv[1:len(sys.argv)])#, shell=False, check=False)
else:
    print("Unknown program to run: "+ sys.argv[1]+"\n")
    print("Admissible program names are: "+str(prog_names) +"\n\n")
