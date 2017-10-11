import sys,script_tree,random,time,math,string

MINIMUM_SUPPORT_WITHIN_A_FAMILY = 0.05
MINIMUM_FAMILY_SIZE = 5
#MAXIMUM_DIFFERENCE_WITH_PROFILE = 6000
OUTPUT_FILE_REC = "constraints_from_recs"
OUTPUT_FILE_TRF = "constraints_from_transfers"

parameters = sys.argv[1:]

species_tree = script_tree.readTree(open(parameters[0],"r").readline())
extant_species = script_tree.getLeavesNames(species_tree)
gene_tree_file = parameters[1]
gene_trees = []
for gt in open(gene_tree_file,"r").readlines():
  gene_trees.append(gt.strip())
nodes = script_tree.getNodes(species_tree)
for n in nodes:
  script_tree.setLength(species_tree,n,1.0)
  
  
  
def distance_from(x,y):
  nodes = script_tree.getNodes(species_tree)
  for n in nodes:
    if script_tree.isLeaf(species_tree,n):
      name = script_tree.getName(species_tree,n)
    else:
      name = script_tree.getBootstrap(species_tree,n)
    if name == x:
      A = n
    if name == y:
      B = n
  return script_tree.distanceFrom(species_tree,A,B)

def parent(x):
  #print x
  nodes = script_tree.getNodes(species_tree)
  for n in nodes:
    if script_tree.isLeaf(species_tree,n):
      name = script_tree.getName(species_tree,n)
    else:
      name = script_tree.getBootstrap(species_tree,n)
    if name == x:
      if script_tree.isRoot(species_tree,n):
	result = -1
      else:
	result = script_tree.getBootstrap(species_tree,script_tree.getParent(species_tree,n))
  return result

def difference_symetrique(A,B):
  result = 0
  for a in A:
    if not a in B:
      result = result + 1
  for a in B:
    if not a in A:
      result = result + 1
  return result
  

def donnor_search(t,n):
  current_node = n
  result = -1
  while result == -1 and not script_tree.isRoot(t,current_node):
    current_node = script_tree.getParent(t,current_node)
    annot = script_tree.getBootstrap(t,current_node)
    event = annot[1:].split(".")[0]
    if event[:2] != "T@" and event[:2] != "D@":
      result = event
    if annot.find("@T") >= 0:
      current_node = script_tree.getRoot(t)
      result = -1
  return result
    
  
def receptor_search(t,n):
  if script_tree.isLeaf(t,n):
    result = []
  else:
    annot = script_tree.getBootstrap(t,n)
    if annot.find("T@") >= 0:
      result = []
    elif annot.split(".")[1][:2] == "D@":
      children = script_tree.getChildren(t,n)
      result = receptor_search(t,children[0]) + receptor_search(t,children[1])
    else:
      result = [annot.split(".")[1]]
  #print result
  return result
    
#output_distr = open("distrib_diff","w")
output_rec = open(OUTPUT_FILE_REC,'w')
output_rec.write("#family,older,younger,support\n")
output_trf = open(OUTPUT_FILE_TRF,'w')
output_trf.write("#family,older,younger,support\n")
for gt in gene_trees:
  file_rec = open(gt,"r").readlines()
  print gt,
  constraints_rec = {}
  constraints_trf = {}
  i = 0
  while i < len(file_rec):
    words = file_rec[i].split()
    #print words,len(words)
    if len(words) > 1 and words[1] == "reconciled":
      number = int(words[0])
      for j in range(i+2,i+2+number):
	#print file_rec[j]
	tree = script_tree.readTree(file_rec[j])  # one reconciliation
	root = script_tree.getRoot(tree)
	#annot = script_tree.getBootstrap(tree,root)
	#event = annot.split(".")[-1]
	#if event[:2] == "T@" or event[:2] == "D@":
	  #origination = event[2:].split("->")[0]
	#else:
	  #origination = event              # detect origination species
	#nodes = script_tree.getNodes(species_tree)
	##for n in nodes:
	  #if script_tree.isLeaf(species_tree,n):
	    #name = script_tree.getName(species_tree,n)
	    #if name == origination:
	      #species_under_origin = [name]
	  #else:
	    #name = script_tree.getBootstrap(species_tree,n)
	    #if name == origination:
	      #leaves = script_tree.getLeaves(species_tree,n)
	      #species_under_origin = []
	      #for l in leaves:
		#species_under_origin.append(script_tree.getName(species_tree,l))
	leaves = script_tree.getLeaves(tree,root)
	#species_present = []
	#for l in leaves:
	  #if script_tree.getName(tree,l).split("_")[0] not in species_present:
	    #species_present.append(script_tree.getName(tree,l).split("_")[0])
	if len(leaves) > MINIMUM_FAMILY_SIZE:# and difference_symetrique(species_present,species_under_origin) < MAXIMUM_DIFFERENCE_WITH_PROFILE:
	  #print script_tree.writeTree(tree,root,False)
	  #print origination,len(leaves),len(species_present),len(species_under_origin),
	  #output_distr.write(str(difference_symetrique(species_present,species_under_origin))+"\n")
	  nodes = script_tree.getNodes(tree)
	  for n in nodes:
	    if script_tree.isLeaf(tree,n):
	      annot = script_tree.getName(tree,n)
	    else:
	      annot = script_tree.getBootstrap(tree,n)
	    #if script_tree.isLeaf(tree,n) and annot.find("T")>=0:
	      #print annot
	    # annot: evenement
	    events = annot.split(".")
	    if events[0] == "":
	      del events[0]
	    for e in events:
	      if e[:2] == "T@":
		# on a touve un transfert
		#print e
		donnor = e[2:].split("->")[0]
		receptor = e[2:].split("->")[1]
		#print annot,donnor,receptor
		# 1/ calculer la contrainte selon les transferts
		if parent(donnor) != -1 and not (receptor in script_tree.getLeavesNames(species_tree)):
		  c = str(parent(donnor))+","+str(receptor)
		  if not constraints_trf.has_key(c):
		    constraints_trf[c] = 0
		  constraints_trf[c] = constraints_trf[c] + 1
		# 2/ calculer la contrainte selon les reconciliations  
		if script_tree.isLeaf(tree,n):
		  receptor = []
		else:
		  if events[0][:2] != "T@":
		    if events[0][:2] == "D@":
		      species = events[0][2:]
		    else:
		      species = events[0]
		    receptor = [species]
		  else:
		    children = script_tree.getChildren(tree,n)
		    child1 = children[0]
		    if script_tree.isLeaf(tree,child1):
		      annot_child1 = script_tree.getName(tree,child1).split(".")[-1].split("_")[0]
		    else:
		      annot_child1 = script_tree.getBootstrap(tree,child1).split(".")[-1]
		    if annot_child1[:2] == "T@":
		      annot_child1 = annot_child1[2:].split("->")[0]
		    if annot_child1[:2] == "D@":
		      annot_child1 = annot_child1[2:]
		    child2 = children[1]
		    if script_tree.isLeaf(tree,child2):
		      annot_child2 = script_tree.getName(tree,child2).split(".")[-1].split("_")[0]
		    else:
		      annot_child2 = script_tree.getBootstrap(tree,child2).split(".")[-1]
		    if annot_child2[:2] == "T@":
		      annot_child2 = annot_child2[2:].split("->")[0]
		    if annot_child2[:2] == "D@":
		      annot_child2 = annot_child2[2:]
		    #print annot,annot_child1,annot_child2
		    if (receptor == annot_child1 and donnor == annot_child2):
		      receptor = receptor_search(tree,child1)
		    elif (donnor == annot_child1 and receptor == annot_child2):
		      receptor = receptor_search(tree,child2)
		    else:
		      receptor = []
		      #print "rien",donnor,receptor
		donnor = donnor_search(tree,n)
		#print donnor,receptor
		#if receptor != []:
		  #print annot,donnor,receptor,len(leaves)    
		for r in receptor:
		  if donnor != -1 and not r in extant_species:
		    if not constraints_rec.has_key(str(donnor)+","+str(r)):
		      constraints_rec[str(donnor)+","+str(r)] = 0
		    constraints_rec[str(donnor)+","+str(r)] = constraints_rec[str(donnor)+","+str(r)] + 1
      i = i + number
    i = i + 1

  for k in constraints_trf.keys():
    donnor_receptor = k
    receptor_donnor = donnor_receptor.split(",")[1]+","+donnor_receptor.split(",")[0]
    if constraints_trf.has_key(receptor_donnor) and constraints_trf[receptor_donnor] <= constraints_trf[donnor_receptor]:
      constraints_trf[donnor_receptor] = constraints_trf[donnor_receptor] - constraints_trf[receptor_donnor]
      constraints_trf[receptor_donnor] = 0

  for k in constraints_rec.keys():
    donnor_receptor = k
    receptor_donnor = donnor_receptor.split(",")[1]+","+donnor_receptor.split(",")[0]
    if constraints_rec.has_key(receptor_donnor) and constraints_rec[receptor_donnor] <= constraints_rec[donnor_receptor]:
      constraints_rec[donnor_receptor] = constraints_rec[donnor_receptor] - constraints_rec[receptor_donnor]
      constraints_rec[receptor_donnor] = 0

  print len(constraints_rec.keys()),sum(constraints_rec.values())
  for k in constraints_rec.keys():
    if constraints_rec[k] / float(number) > MINIMUM_SUPPORT_WITHIN_A_FAMILY:
      output_rec.write(str(gt).split("/")[-1]+","+k+","+str(constraints_rec[k]/float(number))+","+str(distance_from(k.split(",")[1],k.split(",")[0]))+"\n")
  for k in constraints_trf.keys():
    if constraints_trf[k] / float(number) > MINIMUM_SUPPORT_WITHIN_A_FAMILY:
        donnor = k.split(",")[0]
        receptor = k.split(",")[1]
        output_trf.write(str(gt).split("/")[-1]+","+donnor+","+receptor+","+str(constraints_trf[k]/float(number))+","+str(distance_from(k.split(",")[1],k.split(",")[0]))+"\n")
