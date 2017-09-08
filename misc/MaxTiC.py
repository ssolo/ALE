# MaxTiC: Fast ranking of a phylogenetic tree by Maximum Time Consistency with lateral gene transfers
#
#
# python code written by Eric Tannier, Inria
# Using ideas, comments, suggestions from Cedric Chauve, Akbar Rafiey, Adrian A. Davin, Celine Scornavacca, Philippe Veber, Bastien Boussau, Gergely J Szollosi, Vincent Daubin
#
# Software distributed under the cecill licence, rights and permissions are described here:
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#
# Bug reports, suggestions or help request for usage of MaxTic should be sent to eric.tannier@inria.fr
#
# If you use this software please cite
# MaxTiC: Fast ranking of a phylogenetic tree by Maximum Time Consistency with lateral gene transfers, Biorxiv doi.org/10.1101/127548 


import sys,script_tree,random,time,math

MAX_NUMBER = 10000000000
time0 = time.time()

parameters = sys.argv[1:]

if len(parameters) < 2:
  print "usage: python MaxTiC.py species_tree_file constraints_file [ls=LOCAL_SEARCH] [t=TEMPERATURE] [r=RANDOMIZATION_TYPE] [d=MIN_DISTANCE] [rd=RANDOM_TREES] [ts=THRESHOLD_CONSTRAINTS]"
  exit()
  
#RANDOMIZATION_TYPE: 0 = take data as it is, 1 = keep nodes, randomize direction, 2 = randomize nodes

TEMPERATURE = 0.001
RANDOM_TYPE = 0
TIME_FOR_SEARCH = 0
NUMBER_OF_REPEATS = 1
MIN_TRANSFER_DIST = 0
RANDOM = 0
THRESHOLD_CONSTRAINTS = 0.0
for p in parameters[2:]:
  words = p.split("=")
  if words[0] == "t":
    TEMPERATURE = float(words[1])
  elif words[0] == "r":
    RANDOM_TYPE = int(words[1])
  elif words[0] == "ls":
    TIME_FOR_SEARCH = float(words[1])
  elif words[0] == "d":
    MIN_TRANSFER_DIST = int(words[1])
  elif words[0] == "rd":
    RANDOM = int(words[1])
  elif words[0] == "ts":
    THRESHOLD_CONSTRAINTS = float(words[1])
  else:
    print "unused parameter (bad format):",p


name_constraint_file = parameters[1]
tree = script_tree.readTree(open(parameters[0],"r").readline())
constraints = open(name_constraint_file,"r").readlines()
root = script_tree.getRoot(tree)

print name_constraint_file


def maximum_spearman(n):
  if n % 2 == 0:
    m = n /2
    return 2 * m * m
  else:
    m = (n-1)/2
    return 2 * m * m + 2 * m
  
#def maximum_distance(tree,root):
  #result = [[root],[root]]
  #if not script_tree.isLeaf(tree,root):
    #child1 = script_tree.getChildren(tree,root)[0]
    #child2 = script_tree.getChildren(tree,root)[1]
    #orders1 = maximum_distance(tree,child1)
    #orders2 = maximum_distance(tree,child2)
    #result[0] = result[0] + orders1[0]
    #result[0] = result[0] + orders2[1]
    #result[1] = result[1] + orders2[0]
    #result[1] = result[1] + orders1[1]
  #return result

def maximum_distance(tree,root):
  c1 = script_tree.getChildren(tree,root)[0]
  c2 = script_tree.getChildren(tree,root)[1]
  name = script_tree.getBootstrap(tree,root)
  if script_tree.isLeaf(tree,c1) and script_tree.isLeaf(tree,c2):
    return [[name],[name]]
  elif script_tree.isLeaf(tree,c1):
    oo = maximum_distance(tree,c2)
    return [[name] + oo[0],[name] + oo[1]]
  elif script_tree.isLeaf(tree,c2):
    oo = maximum_distance(tree,c1)
    return [[name] + oo[0],[name] + oo[1]]
  else:
    orders1 = maximum_distance(tree,c1)
    orders2 = maximum_distance(tree,c2)
    return [[name] + orders1[0] + orders2[0],[name] + orders2[1] + orders1[1]]

def spearman_distance(A,B):
  result = 0.0
  Binv = {}
  for i in range(len(B)):
    Binv[B[i]] = i
  #print len(Binv.keys()),len(A)
  for i in range(len(A)):
    j = Binv[A[i]]
    #print i,j,A[i]
    result = result + abs(i-j)
  #print maximum_spearman(len(A)),result
  return result

def spearman_similarity(A,B):
  m = maximum_distance(tree,script_tree.getRoot(tree))
  #print m,A,B
  #print spearman_distance(m[0],m[1]),spearman_distance(A,B)
  #print (spearman_distance(m[0],m[1]) - spearman_distance(A,B))/(spearman_distance(m[0],m[1]))
  return (spearman_distance(m[0],m[1]) - spearman_distance(A,B))/(spearman_distance(m[0],m[1]))

def kendall_distance(A,B):
  result = 0.0
  Ainv = {}
  for i in range(len(A)):
    Ainv[A[i]] = i
  Binv = {}
  for i in range(len(B)):
    Binv[B[i]] = i
  for i in range(len(A)):
    for j in range(i+1,len(A)):
      if Binv[A[i]] > Binv[A[j]]:
	result = result + 1
  return result

def kendall_similarity(A,B):
  m = maximum_distance(tree,script_tree.getRoot(tree))
  #print m
  #print kendall_distance(m[0],m[1]) , kendall_distance(A,B)
  return (kendall_distance(m[0],m[1]) - kendall_distance(A,B))/(kendall_distance(m[0],m[1]))

def similarity(A,B):
  #return spearman_similarity(A,B)
  return kendall_similarity(A,B)


def path(g,a,b):  # is there a directed path in graph g from a to b, true if a=b
        marques = [a]
        pile = [a]
        while len(pile) > 0 and (not (b in marques)):
                sommet = pile[-1]
                del pile[-1]
                voisins = g[sommet]
                for v in voisins:
                        if not v in marques:
                                marques.append(v)
                                pile.append(v)
        if b in marques:
                return True
        else:
                return False



def value(order):
  index = {}
  for i in range(len(order)):
    index[order[i]] = i
  result = 0
  for e in edge_keys:
    sommets = e.split(",")
    if index.has_key(sommets[0]) and index.has_key(sommets[1]) and index[sommets[0]] > index[sommets[1]]:
      result = result + edge[e]
  return result

def tree_from_order(order):
  nodes = script_tree.getNodes(tree)
  for n in nodes:
    if not script_tree.isRoot(tree,n):
      if script_tree.isLeaf(tree,n):
	index = len(order)
      else:
	index = order.index(script_tree.getBootstrap(tree,n))
      index_parent = order.index(script_tree.getBootstrap(tree,script_tree.getParent(tree,n)))
      script_tree.setLength(tree,n,index-index_parent)

  return script_tree.writeTree(tree,root,False)

def order_from_graph(graph):
  result = []
  marques = {}
  for l in script_tree.getLeavesNames(tree):
    marques[l] = 0
  while len(result) < len(graph.keys()) - len(script_tree.getLeavesNames(tree)):
    keys = graph.keys()
    i = 0
    while marques.has_key(keys[i]):
      i = i + 1
    current = keys[i]
    suivant = True
    while suivant:
      suivant = False
      for v in graph[current]:
	if not marques.has_key(v):
	  suivant = True
	  current = v
    result.append(current)
    marques[current] = 0
  result.reverse()
  return result
  
  while len(graph[current]) > 0:
    current = graph[current][0]
    result.append(current)
  return result

def optimisation_locale(order,duration):
  #sortie_tree = open("MT_output_tree_sample","w")
  #sortie_val = open("MT_output_score_sample","w")
  time0=time.time()
  current = value(order)
  best = current
  best_order = list(order)
  sample = []
  #print ref
  while time.time() - time0 < duration:
    i = int(random.random()*(len(order)))
    j = int(random.random()*(len(order)))
    if i != j:
      essai = list(order)
      a = min(i,j)
      b = max(i,j)
      temp = essai[a]
      essai[a:b] = essai[a+1:b+1]
      essai[b] = temp
      v = value(essai)
      if v < MAX_NUMBER:
	#sortie_tree.write(tree_from_order(order)+" "+str(current)+" "+str(TEMPERATURE)+"\n")
	#sortie_val.write(str(current)+"\n")
	if v <= current:
	  metropolis_ratio = 1
	  #TEMPERATURE = min(0.01,TEMPERATURE/10)
	else:
	  metropolis_ratio = math.exp((current-v)/TEMPERATURE)
	coin = random.random()
	if coin < metropolis_ratio:
	  order = essai
	  current = v
	if v < best:
          print "better solution",v  
	  best = v
	  best_order = list(order)
	  sample = []
	if v == best and len(sample) < 10000:
	  sample.append(similarity(order_input,order))
	  
  #sortie_val = open("MT_output_score_sample","w")
  #for s in sample:
    #sortie_val.write(str(s)+"\n")
  return best_order
    

def edgeweights(element,elements):
  result = 0
  for s in degre_entrant[element]:
    if s in elements:
      if edge.has_key(s+","+element):
	result = result + edge[s+","+element]
  return result

def mix(order1,order2):
  order = []
  cost = [[0]*(len(order2)+1)] 
  for i in range(1,len(order1)+1):
    cost.append([0]*(len(order2)+1))
  back = [["j"]*(len(order2)+1)] 
  for i in range(1,len(order1)+1):
    back.append(["i"]+[""]*(len(order2)))
  for i in range(1,len(order1)+1):
    for j in range(1,len(order2)+1):
      #print i,j,len(order1),len(order2)
      value1 = cost[i-1][j] + edgeweights(order1[i-1],order2[0:j])
      value2 = cost[i][j-1] + edgeweights(order2[j-1],order1[0:i])
      if value1 == value2:
        x = random.random()
        if x < 0.5:
          cost[i][j] = value1
          back[i][j] = "i"
        else:
          cost[i][j] = value2
          back[i][j] = "j"
      elif value1 < value2:
        cost[i][j] = value1
        back[i][j] = "i"
      else:
        cost[i][j] = value2
        back[i][j] = "j"
  i = len(order1)
  j = len(order2)
  while i>0 or j>0:
    if back[i][j] == "i":
      #print i,len(order1)
      order.append(order1[i-1])
      i = i - 1
    else:
      order.append(order2[j-1])
      j = j - 1
  order.reverse()
  return order
     
def opt(tree,root):
  if (script_tree.isLeaf(tree,script_tree.getChildren(tree,root)[0]) and
      script_tree.isLeaf(tree,script_tree.getChildren(tree,root)[1])):
    return [script_tree.getBootstrap(tree,root)]
  elif script_tree.isLeaf(tree,script_tree.getChildren(tree,root)[0]):
    return opt(tree,script_tree.getChildren(tree,root)[1])+[script_tree.getBootstrap(tree,root)]
  elif script_tree.isLeaf(tree,script_tree.getChildren(tree,root)[1]):
    return opt(tree,script_tree.getChildren(tree,root)[0])+[script_tree.getBootstrap(tree,root)]
  else:
    order1 = opt(tree,script_tree.getChildren(tree,root)[0])
    order2 = opt(tree,script_tree.getChildren(tree,root)[1])
    #print "level",len(order1)+len(order2),order1,order2
    return mix(order1,order2) + [script_tree.getBootstrap(tree,root)]
  

def random_order(tree,root):
  c1 = script_tree.getChildren(tree,root)[0]
  c2 = script_tree.getChildren(tree,root)[1]
  name = script_tree.getBootstrap(tree,root)
  if script_tree.isLeaf(tree,c1) and script_tree.isLeaf(tree,c2):
    return [name]
  elif script_tree.isLeaf(tree,c1):
    return [name] + random_order(tree,c2)
  elif script_tree.isLeaf(tree,c2):
    return [name] + random_order(tree,c1)
  else:
    order1 = random_order(tree,c1)
    order2 = random_order(tree,c2)
    #if script_tree.isRoot(tree,root):
      #print order1,order2
    pos = range(len(order1)+len(order2))
    for i in range(len(order2)):
      index = int(random.random()*(len(pos)))
      del pos[index]
    #if script_tree.isRoot(tree,root):
      #print pos
    order = [name]
    previous = 0
    for i in pos:
      for j in range(i - previous):
	#if script_tree.isRoot(tree,root):
	  #print order2,i,previous
	order.append(order2[0])
	del order2[0]
      order.append(order1[0])
      del order1[0]
      previous = i+1
	
    #if script_tree.isRoot(tree,root):
      #print order1,order2,order
    return order + order2


def order_from_tree(tree):
  order = []
  nodes = script_tree.getNodes(tree)
  root = script_tree.getRoot(tree)
  for n in nodes:
    if not script_tree.isLeaf(tree,n):
      order.append(n)
  #print order
  order.sort(lambda x,y: cmp(script_tree.distanceFrom(tree,x,root),script_tree.distanceFrom(tree,y,root)))
  for i in range(len(order)):
    order[i] = script_tree.getBootstrap(tree,order[i])
  return order

# useful variables
degre_entrant = {}
edge = {}
graph = {}
nodes = script_tree.getNodes(tree)
leaves = script_tree.getLeavesNames(tree)
root = script_tree.getRoot(tree)
internal_nodes = []

# initialise from nodes
for n in nodes:
  if (not script_tree.isLeaf(tree,n)):
    node = script_tree.getBootstrap(tree,n)
    graph[node] = []
    degre_entrant[node] = []
    internal_nodes.append(node)
  else:
    node = script_tree.getName(tree,n)
    graph[node] = []

# initialise from branches
for n in nodes:
  if (not script_tree.isLeaf(tree,n)):
    if (not script_tree.isRoot(tree,n)):
      parent = script_tree.getBootstrap(tree,script_tree.getParent(tree,n))
      child = script_tree.getBootstrap(tree,n)
      graph[parent].append(child)
      edge[parent+","+child] = MAX_NUMBER
  else:
    parent = script_tree.getBootstrap(tree,script_tree.getParent(tree,n))
    child = script_tree.getName(tree,n)
    graph[parent].append(child)

print "tree with ",len(internal_nodes),"internal nodes"


# parse transfers
#by_family = {}
for line in constraints:
  #print line
  if (len(line.split()) >= 2 or len(line.split(",")) >= 2) and line.find("FRQ") <0 and line[0] != "#":
    #print line
    # compute first and second, the donnor and receptor
    if line.find(",")>=0:
      words = line.strip().split(",")
      words = words[1:]
    else:
      words = line.split()
    if len(words) <= 3 or float(words[3]) > MIN_TRANSFER_DIST:
      if len(words) > 2:
	weight = float(words[2])
      else:
	weight = 1.0
      if RANDOM_TYPE == 0:
	first = words[0]
	second = words[1]
      elif RANDOM_TYPE == 1:
	first = words[0]
	second = words[1]
	if random.random() < 0.5:
	  tmp = first
	  first = second
	  second = tmp
      else:
	first = internal_nodes[int(random.random()*len(internal_nodes))]
	second=first
	while second == first:
	  second=internal_nodes[int(random.random()*len(internal_nodes))]
	#print first,second
      if first == "None":
	first = script_tree.getBootstrap(tree,root)
      if second == "None":
	second = script_tree.getBootstrap(tree,root)
      # first and second have been chosen or read
      if True or (first in internal_nodes and second in internal_nodes):
	key = first+","+second
	if not edge.has_key(key):
            edge[key] = 0
        if not edge[key] >= MAX_NUMBER:
            edge[key] = edge[key] + weight
            #if by_family.has_key(family):
            #by_family[family] = by_family[family] + weight
  else:
    print line.strip(),"                       ...ignored"
    
#sortie = open("distribution_by_family","w")
#keys = by_family.keys()
#keys.sort()
#for f in keys:
  #sortie.write(f+" "+str(by_family[f])+"\n")


edge_keys = edge.keys()
edge_keys.sort(lambda x,y: cmp(edge[x],edge[y]))


total_transfers = 0.0
for k in edge_keys:
    if edge[k] < MAX_NUMBER:
        total_transfers = total_transfers + edge[k]

sub_total = 0.0
while sub_total < total_transfers*THRESHOLD_CONSTRAINTS:
    sub_total = sub_total + edge[edge_keys[0]]
    del edge[edge_keys[0]]
    del edge_keys[0]

# remove under THRESHOLD_CONSTRAINTS
total_transfers = 0.0
for k in edge.keys():
    if edge[k] < MAX_NUMBER:
        first = k.split(",")[0]
        second = k.split(",")[1]
        if degre_entrant.has_key(second):
            degre_entrant[second].append(first)
        total_transfers = total_transfers + edge[k]

print total_transfers,"total weight of constraints from transfers"

# remove uninformative
trivial_conflict = 0
to_itself = 0
to_leaf = 0
from_leaf = 0
to_desc = 0
to_anc = 0
uninformative = 0
edge_keys = edge.keys()
edge_keys.sort(lambda x,y: cmp(edge[y],edge[x]))
for e in edge_keys:
  first = e.split(",")[0]
  second = e.split(",")[1]
  if (first == second) and edge[e] < MAX_NUMBER:
    uninformative = uninformative + edge[e]
    to_itself = to_itself + edge[e]
  elif path(graph,first,second) and edge[e] < MAX_NUMBER:
    uninformative = uninformative + edge[e]
    to_desc = to_desc + edge[e]
  elif (second in leaves) and edge[e] < MAX_NUMBER:
    uninformative = uninformative + edge[e]
    to_leaf = to_leaf + edge[e]
  elif (first in leaves) and edge[e] < MAX_NUMBER:
    uninformative = uninformative + edge[e]
    from_leaf = from_leaf + edge[e]
    #print first,second
  elif path(graph,second,first):
     uninformative = uninformative + edge[e]
     trivial_conflict = trivial_conflict + edge[e]
     #print first,second
     to_anc = to_anc + edge[e]

# output informative
edge_keys = edge.keys()
edge_keys.sort(lambda x,y: cmp(edge[y],edge[x]))


# substract opposite
for e in edge_keys:
  words = e.split(",")
  opposite = words[1]+","+words[0]
  #print e,opposite
  if edge[e] < MAX_NUMBER and edge.has_key(opposite) and edge[opposite] < MAX_NUMBER and edge[e]>= edge[opposite]:
    #edge[e] = edge[e] - edge[opposite]
    #total_transfers = total_transfers
    trivial_conflict = trivial_conflict + edge[opposite]
    #print e
    #edge[opposite] = 0

for e in edge.keys():
  if edge[e] == 0:
    del edge[e]
   
# output informative
edge_keys = edge.keys()
edge_keys.sort(lambda x,y: cmp(edge[y],edge[x]))
sortie = open(name_constraint_file+"_MT_output_filtered_list_of_weighted_informative_constraints","w")
for k in edge_keys:
  if edge[k] < MAX_NUMBER:
    sortie.write(k+" "+str(edge[k])+"\n")
sortie.close()

print uninformative,"uninformative (",to_desc,"to a descendant,",to_leaf,"to a leaf",to_anc,"to an ancestor",to_itself,"to itself",
print str(int(uninformative*100/(total_transfers+uninformative)))+"%)"
if from_leaf > 0:
    print "WARNING: there are constraints from leaves, which cannot be met since this program only orders internal nodes, they are ignored in the scores"
print trivial_conflict,"trivially conflicting constraints ("+str(int(trivial_conflict*100/(total_transfers)))+"%) (descendant to ancestor or trivial cycle)"
#print total_transfers,"remaining informative constraints from transfers (removing trivial cycles)"

order_input = order_from_tree(tree)
value_input = value(order_input)
print "value of the order given by the input tree",value_input,"("+str((value_input*100/total_transfers))+"%)"

#print "("+str(trivial_conflict+value(order_input))+" or "+str(int((trivial_conflict+value(order_input))*100/(total_transfers+2*trivial_conflict)))+"% of total informative constraints)"

#print order_input

edge_keys.sort(lambda x,y: cmp(edge[y],edge[x]))
rejected = 0
for e in edge_keys:
  #print e
  words = e.split(",")
  if words[0] != words[1] and path(graph,words[1],words[0]):
    rejected = rejected + edge[e]
  else:
    graph[words[0]].append(words[1])
#print "greedy solution",rejected,"("+str(int(rejected*100/total_transfers))+"%)"
#print "("+str(int((trivial_conflict+rejected)*100/(total_transfers+trivial_conflict)))+"% of total informative constraints)"

order_greedy = order_from_graph(graph)
value_greedy = value(order_greedy)

print "value of the greedy heuristic",value_greedy,"("+str((value_greedy*100/total_transfers))+"%)"

edge_keys = edge.keys()
edge_keys.sort(lambda x,y: cmp(edge[y],edge[x]))

time0 = time.time()
order_heuristic = opt(tree,root)
order_heuristic.reverse()
value_heuristic = value(order_heuristic)
print "value of the mixing heuristic:",value_heuristic,"("+str((value_heuristic*100/total_transfers))+"%)"
#print "("+str(trivial_conflict+value(order))+" or ",str(int((trivial_conflict+value(order))*100/(total_transfers+2*trivial_conflict)))+"% of total informative constraints) in",time.time()-time0,"secs"

#if value_input <= value_greedy and value_input <= value_heuristic:
#  order = order_input
#  print "best order is the input tree"
#if value_greedy <= value_heuristic and value_greedy <= value_input:
  #order = order_greedy
  #print "best order is the greedy heuristic"
#if value_heuristic <= value_greedy and value_heuristic <= value_input:
  #order = order_heuristic
  #print "best order is the mixing heuristic"
if value_greedy <= value_heuristic:
    order = order_greedy
    print "best order is the greedy heuristic"
else:
    order = order_heuristic
    print "best order is the mixing heuristic"
print tree_from_order(order)
#print order
#print order_input
print "Similarity of the best order compared with the input order", similarity(order_input,order)


if RANDOM>0:
  sortie = open("distribution_random","w")
  minimum = MAX_NUMBER
  maximum = 0
  pvalue = 0
  kendall_tau_min = MAX_NUMBER
  kendall_tau_max = 0
  kendall_tau_pvalue = 0
  for i in range(RANDOM):
    rorder = random_order(tree,root)
    val = value(rorder)
    sortie.write(str(val)+" ")
    minimum = min(val,minimum)
    maximum = max(val,maximum)
    if val <= value(order):
      pvalue = pvalue + 1
    val = similarity(rorder,order_input)
    sortie.write(str(val) + "\n")
    kendall_tau_min = min(val,kendall_tau_min)
    kendall_tau_max = max(val,kendall_tau_max)
    if val >= similarity(order_input,order):
      kendall_tau_pvalue = kendall_tau_pvalue + 1
    total=RANDOM
    position=i+1
    sys.stdout.write("\r generating "+str(RANDOM)+" ranked trees ["+
                         "="*int(float(position)*30/float(total))+
                         " "*(30-int(float(position)*30/float(total)))+
                         "]  "+str(int(float(position)*100/float(total)))+
                         "% executed, after "+str(int(time.time()-time0))+" sec")
  print
  print "values from ",RANDOM," random orders",minimum,maximum
  print "pvalue of the found order:",float(pvalue)/RANDOM
  print "similarity values from ",RANDOM,"random orders",kendall_tau_min,kendall_tau_max
  print "pvalue of the similarity with the input order:",float(kendall_tau_pvalue)/RANDOM
  


if TIME_FOR_SEARCH > 0:
  print "attempting a local search from the best found order, please wait",TIME_FOR_SEARCH,"seconds"
  order = optimisation_locale(order,TIME_FOR_SEARCH)
  print "after local search", value(order),"rejected"
  print "best found solution",value(order),"("+str((value(order)*100/total_transfers))+"%)"
  print tree_from_order(order)
  print order

sortie = open(name_constraint_file+"_MT_output_list_of_constraints_conflicting_with_best_order","w")
index = {}
for i in range(len(order)):
  index[order[i]] = i
for e in edge_keys:
  sommets = e.split(",")
  if index.has_key(sommets[0]) and index.has_key(sommets[1]) and index[sommets[0]] > index[sommets[1]]:
      sortie.write(e+" "+str(edge[e])+"\n")
    
conflict_with_input = 0.0
total = 0.0
sortie = open(name_constraint_file+"_MT_output_partial_order","w")
for e in edge_keys:
  sommets = e.split(",")
  if index.has_key(sommets[0]) and index.has_key(sommets[1]) and index[sommets[0]] < index[sommets[1]] and edge[e] < 100000:
    total = total + edge[e]
    if order_input.index(sommets[0]) < order_input.index(sommets[1]):
      sortie.write(str(sommets[0])+" "+str(sommets[1])+" "+str(edge[e])+" black\n")
    else:
      sortie.write(str(sommets[0])+" "+str(sommets[1])+" "+str(edge[e])+" green\n")
      conflict_with_input = conflict_with_input + edge[e]

print conflict_with_input/total," constraints in agreement with the best tree in conflict with input tree"

print "Similarity of the order compared with the input order", similarity(order_input,order)

#orders = max_dist_orders(tree,root)
#print orders
#print "max dist orders",score_order(orders[0],orders[1])
