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
 #   p = subprocess.Popen(sys.argv[1:len(sys.argv)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
 #   (out, err) = p.communicate()
 #   if err:
 #       print(err)
 #   print(out)
#self.logger.error(err)
#        return (out, err)
    subprocess.call(sys.argv[1:len(sys.argv)])#, shell=False, check=False)
#    subprocess.call("ls -ltrh", shell=True)#, shell=False, check=False)
#    print(subprocess.check_output(["ls", "-ltrh"]))
#    print(subprocess.check_output(["touch", "toto"]))
#    print(subprocess.check_output(["ls", "-ltrh"]))
else:
    print("Unknown program to run: "+ sys.argv[1]+"\n")
    print("Admissible program names are: "+str(prog_names) +"\n\n")

print ("Thank you for using ALEsuite.")
