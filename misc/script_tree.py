###############################################################
####### PROCESSING OF TREES ###################################
###############################################################

# structure of the tree:
# 0: name, 1: parent, 2: tab of children, 3: length, 4: isdup, 5:species, 6:bootstrap , 7: bppnumber, 8: ND

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
  
def getbppnumber(tree,node):
	return tree[node][7]

def writeBootstrap(tree,node,value):
	tree[node][6] = value

def getBootstrap(tree,node):
	return tree[node][6]

def addNode(tree):
    id_node = 0
    while tree.has_key(id_node):
        id_node = id_node + 1
    tree[id_node] = ["N"+str(id_node),-1,[],0,"","",""]
    return id_node

def getAncestor(tree):
    if tree.has_key("ancestor"):
        return tree["ancestor"]
    else:
        return -1

def setAncestor(tree,node):
    tree["ancestor"] = node
    
def getLength(tree,node):
    return tree[node][3]

def setLength(tree,node,l):
    tree[node][3] = l

def getName(tree,node):
    return tree[node][0]

def setName(tree,node,name):
    tree[node][0] = name

def getSpecies(tree,node):
    return tree[node][5]

def writeSpecies(tree,node,annot):
    tree[node][5] = annot

def getNodes(tree):
    clefs = tree.keys()
    c = 0
    while c < len(clefs):
        if (clefs[c] == "sequence" or
            clefs[c] == "ancestor" or
            len(tree[clefs[c]]) == 0):
            del clefs[c]
        else:
            c = c + 1
    return clefs

def getParent(tree,node):
    return tree[node][1]

def getChildNumber(tree,n,c):
  children = getChildren(tree,n)
  if children[0] == c:
    return 0
  else:
    return 1

def setParent(tree,node,p):
  tree[node][1] = p

def addChild(tree,pere,node):
    tree[pere][2].append(node)
    
def removeLeaf(tree,l):
  #root = getRoot(tree)  
  #print "remove",l,writeTree(tree,root,False)
  if isRoot(tree,l):
      del tree[l]
  else:
    pere = getParent(tree,l)
    if isRoot(tree,pere) and len(getChildren(tree,pere)) == 2:
        #print "son of the root"
        b = getBrother(tree,l)
        tree[b][1] = -1
        del tree[pere]
        del tree[l]
    elif len(getChildren(tree,pere)) == 2:
        b = getBrother(tree,l)
        grandpere = getParent(tree,pere)
        setParent(tree,b,grandpere)
        number = getChildNumber(tree,grandpere,pere)
        setChild(tree,grandpere,number,b)
        tree[b][3] = tree[b][3]+tree[pere][3]
        del tree[pere]
        del tree[l]
    elif isRoot(tree,pere) and len(getChildren(tree,pere)) == 1:
        del tree[l]
        del tree[pere]
    elif len(getChildren(tree,pere)) > 2:
        number = getChildNumber(tree,pere,l)
        del tree[pere][2][number]
        del tree[l]

    

def removeNodeAndChildren(tree,node):
    children = list(getChildren(tree,node))
    for child in children:
        removeNodeAndChildren(tree,child)
    removeNode(tree,node)

def removeNode(tree,node):
#    print "effacement du noeud",node
    del tree[node]

def removeChildAndChildren(tree,pere,node):
    numero = 0
    while node != getChild(tree,pere,numero):
        numero = numero + 1
    del tree[pere][2][numero]
    removeNodeAndChildren(tree,node)
    
def removeChild(tree,pere,node):
    numero = 0
    while node != getChild(tree,pere,numero):
        numero = numero + 1
    del tree[pere][2][numero]  
    removeNode(tree,node)

def getChild(tree,node,k):
	return tree[node][2][k]

def setChild(tree,node,k,c):
	tree[node][2][k] = c

def getNumberOfChildren(tree,node):
    return len(tree[node][2])

def getChildren(tree,node):
    return tree[node][2]

def getBrother(tree,node):
    anc = getParent(tree,node)
    if (getChild(tree,anc,0) == node):
        return getChild(tree,anc,1)
    else:
        return getChild(tree,anc,0)

def isLeaf(tree,node):
    return (len(getChildren(tree,node)) == 0)

def isRoot(tree,node):
    return (tree[node][1] == -1)

def isDup(tree,node):
    return (tree[node][4] == "D")
  
def getND(tree,node):
  if tree[node].has_key("ND"):
    return tree[node]["ND"]
  else:
    return ""

def lastCommonAncestor(tree,a,b):
    ancestor = -1
    ancestorsa = [a]
    while not isRoot(tree,a):
        a = getParent(tree,a)
        ancestorsa.append(a)
    ancestorsb = [b]
    while not isRoot(tree,b):
        b = getParent(tree,b)
        ancestorsb.append(b)
#    print ancestorsa,ancestorsb
    while len(ancestorsa) > 0 and len(ancestorsb) > 0 and ancestorsa[-1] == ancestorsb[-1]:
        ancestor = ancestorsa[-1]
        del ancestorsa[-1]
        del ancestorsb[-1]
#    print "ancestor",ancestor
    return ancestor

def distanceFrom(tree,a,b):
    ancestor = lastCommonAncestor(tree,a,b)
    distance = 0
    while a != ancestor:
        #print tree[a]
        distance = distance + tree[a][3]
        a = getParent(tree,a)
    while b != ancestor:
        #print tree[b]
        distance = distance + tree[b][3]
        b = getParent(tree,b)        
    return distance

def getLeaves(tree,a):
#    print "getleaves",a
    if isLeaf(tree,a):
	return [a]
    else:
        #print "non feuille",child1(a),child2(a)
        result = []
        children = list(getChildren(tree,a))
        for child in children:
            result = result + getLeaves(tree,child)
        return result
    
    
def writeTree(tree,a,NHX):
#    print a,tree[a]
    if isLeaf(tree,a):
        if isRoot(tree,a):
            chaine = "("
        else:
            chaine = ""
	chaine = chaine + tree[a][0]
	if tree[a][3] != -1:
		#~ print tree[a][3]
		chaine = chaine + ":" + str(tree[a][3])
        if NHX and tree[a][5] != "":
            chaine = chaine + "[&&NHX:S="+tree[a][5]+"]"
        if isRoot(tree,a):
            chaine = chaine + ")" + str(getBootstrap(tree,a))        
    else:
        chaine = "("
        children = list(getChildren(tree,a))
        for child in children:
            chaine = chaine + writeTree(tree,child,NHX)+","
        chaine = chaine[:-1]+")"+str(getBootstrap(tree,a))
        if (not isRoot(tree,a)) and tree[a][3] != -1:
            chaine = chaine + ":" + str(tree[a][3])  
        if NHX and (tree[a][4] != "" or tree[a][5] != ""):
            chaine = chaine + "[&&NHX:"
            if tree[a][5] != "":
                chaine = chaine + "S="+tree[a][5]
            if tree[a][4] == "D" or  tree[a][4] == "WGD":
                chaine = chaine+":D=Y"
            chaine = chaine + "]"
    if isRoot(tree,a):
        chaine = chaine + ";"
    return chaine

def writeTreeNexus(tree,a,tab):
#    print a,tree[a]
    if isLeaf(tree,a):
        if isRoot(tree,a):
            chaine = "("
        else:
            chaine = ""
	chaine = chaine + tree[a][0]
	chaine = chaine + "[&!color=#"+tab[a]+"]"
	if tree[a][3] != -1:
		#~ print tree[a][3]
		chaine = chaine + ":" + str(tree[a][3])
        if isRoot(tree,a):
            chaine = chaine + ")"
    else:
        chaine = "("
        children = list(getChildren(tree,a))
        for child in children:
            chaine = chaine + writeTreeNexus(tree,child,tab)+","
        chaine = chaine[:-1]+")"
	chaine = chaine + "[&!color=#"+tab[a]+"]"
        if (not isRoot(tree,a)) and tree[a][3] != -1:
            chaine = chaine + ":" + str(tree[a][3])  
    if isRoot(tree,a):
        chaine = chaine + ";"
    return chaine

def getRoot(tree):
    keys = getNodes(tree)
    #print tree
    #print keys
    start = keys[0]
    while (not isRoot(tree,start)):
        start = getParent(tree,start)
    return start

def getNodesBetween(tree,a,b):
    chemin = []
    ancestor = -1
    ancestorsa = []
    while not isRoot(tree,a):
        a = getParent(tree,a)
        ancestorsa.append(a)
    ancestorsb = []
    while not isRoot(tree,b):
        b = getParent(tree,b)
        ancestorsb.append(b)
    while len(ancestorsa) > 0 and len(ancestorsb) > 0 and ancestorsa[-1] == ancestorsb[-1]:
        ancestor = ancestorsa[-1]
        del ancestorsa[-1]
        del ancestorsb[-1]
#    print "ancestor",ancestor
    return ancestorsa+[ancestor]+ancestorsb

def isAncestor(tree,a,b):
        if isRoot(tree,a):
                result = True
        else:
                result = False
        current = b
        while ((not result) and (not isRoot(tree,current))):
                if current == a:
                        result = True
                else:
                        current = getParent(tree,current)
        return result
 
    
def treeCopy(tree):
	result = {}
	for k in tree.keys():
		if k == "ancestor" or k == "sequence":
			result[k] = tree[k]
		else:
			result[k] = [tree[k][0],tree[k][1],list(tree[k][2]),tree[k][3],tree[k][4],tree[k][5],tree[k][6]]
	return result
	

def changeRoot(tree,newRoot): # the new root is between newRoot and its parent NEVER TESTED
	#~ print "changeroot",newRoot,getRoot(tree),getParent(tree,newRoot)
	if (not isRoot(tree,newRoot)) and (not isRoot(tree,getParent(tree,newRoot))):
		#~ print "changeroot"
		root = getRoot(tree)
		new_id = addNode(newtree)
		tree[new_id][2] = [newRoot,getParent(tree,newRoot)]
		tree[newRoot][1] = new_id
		tree[newRoot][3] = tree[newRoot][3]/2
		current = getParent(tree,newRoot)
		prec = new_id
		current_length = tree[newRoot][3]/2
		while getParent(tree,current) != root:
			if current[2][0] == prec:
				tree[current][2][0] = getParent(tree,current)
			else:
				tree[current][2][1] = getParent(tree,current)
			tree[current][1] = prec
			temp = current_length
			current_length = tree[current][3]
			tree[current][3] = temp
			prec = current
			current = getParent(tree,current)
		if current[2][0] == prec:
			tree[current][2][0] = getBrother(tree,current)
		else:
			tree[current][2][1] = getBrother(tree,current)
		tree[current][1] = prec
		temp = current_length
		current_length = tree[current][3]
		tree[current][3] = temp
		tree[getBrother(tree,current)][1] = current
		tree[getBrother(tree,current)][3] = tree[getBrother(tree,current)][3] + current_length
		del tree[root]
			
def SPR(tree,a,b):
	
	#~ print a,b,getName(tree,a),getName(tree,b)
	
	#~ print writeTree(tree,getParent(tree,a),False)
	
	parent = getParent(tree,a)
	great_parent = getParent(tree,getParent(tree,a))
	brother = getBrother(tree,a)
	tree[brother][1] = great_parent
	child = getChildren(tree,great_parent)[0]
	if child == getParent(tree,a):
		tree[great_parent][2][0] = brother
	else:
		tree[great_parent][2][1] = brother
	del tree[parent]
	
	#~ print writeTree(tree,great_parent,False)
	
	parent = getParent(tree,b)
	new_node = addNode(tree)
	tree[new_node][1] = parent
	tree[new_node][2] = [a,b]
	tree[a][1] = new_node
	tree[b][1] = new_node
	child = getChildren(tree,parent)[0]
	if child == b:
		tree[parent][2][0] = new_node
	else:
		tree[parent][2][1] = new_node

	#~ print writeTree(tree,parent,False)

def NNI(tree,node):
  if (not isRoot(tree,node)) and (not isLeaf(tree,node)):
    parent = getParent(tree,node)
    if isRoot(tree,parent):
      brother = getBrother(tree,node)
      if not isLeaf(tree,brother):
	son1 = getChildren(tree,node)[0]
	son2 = getChildren(tree,node)[1]
	son3 = getChildren(tree,brother)[0]
	son4 = getChildren(tree,brother)[1]
	setChild(tree,node,1,son4)
	setChild(tree,brother,1,son2)
	setParent(tree,son2,brother)
	setParent(tree,son4,node)
    else:
      brother = getBrother(tree,node)
      if getChildren(tree,parent)[0] == brother:
	no_brother = 0
      else:
	no_brother = 1
      son1 = getChildren(tree,node)[0]
      setChild(tree,node,0,brother)
      setChild(tree,parent,no_brother,son1)
      setParent(tree,son1,parent)
      setParent(tree,brother,node)
	
	
	


def getLeavesNames(tree):
	result = []
	root = getRoot(tree)
	leaves = getLeaves(tree,root)
	for l in leaves:
		result.append(getName(tree,l))
	return result
	
def unroot(tree):
	nodes = getNodes(tree)
	if len(nodes) > 3:
		root = getRoot(tree)
		children = getChildren(tree,root)
		if len(children) == 2:
			new_root = children[0]
			tree[new_root][1] = -1
			tree[new_root][2].append(children[1])
			tree[children[1]][1] = new_root
			tree[children[1]][3] = tree[new_root][3] + tree[children[1]][3]
			tree[children[1]][6] = max(tree[new_root][6],tree[children[1]][6])
			tree[new_root][3] = -1
			del tree[root]

def contractunsupported(tree,threshold):
        result = 0
	unroot(tree)
	nodes = getNodes(tree)
	#print "begin",len(nodes)
	for n in nodes:
		if isfloat(tree[n][6]):
			tree[n][6] = float(tree[n][6])
		else:
			tree[n][6] = 0.0
		if (not isRoot(tree,n)) and (not isLeaf(tree,n)) and (tree[n][6] < threshold):
			#~ print "CONTRACTION",float(tree[n][6]),threshold,
			parent = getParent(tree,n)
			children = getChildren(tree,n)
			for c in children:
				tree[parent][2].append(c)
				tree[c][1] = parent
			removeChild(tree,parent,n)
			result = result + 1
	#nodes = getNodes(tree)
	#print "end",len(nodes)
	return result
			

def ultrametricize(tree):
	root = getRoot(tree)
	leaves = getLeaves(tree,root)
	maximum = 0
	index = -1
	for l in leaves:
		d = distanceFrom(tree,root,l)
		if d > maximum:
			maximum = d
			index = l
	#~ print getName(tree,l),"maximum",maximum
	i = index
	marque = []
	while i != root:
		marque.append(i)
		i = getParent(tree,i)
	marque.append(root)
	while len(marque) < len(getNodes(tree)):
		#~ print len(marque),len(getNodes(tree))
		maximum_non_marque = 0
		index = -1
		for l in leaves:
			d = distanceFrom(tree,root,l)
			if (d > maximum_non_marque) and (not l in marque):
				maximum_non_marque = d
				index = l
			#~ print getName(tree,l),"distance",distanceFrom(tree,root,l)
		i = index
		distance_from_marque = 0
		while not i in marque:
			distance_from_marque = distance_from_marque + getLength(tree,i)
			i = getParent(tree,i)
		ratio = (maximum - distanceFrom(tree,i,root)) / distance_from_marque
		i = index
		while not i in marque:
			marque.append(i)
			setLength(tree,i,getLength(tree,i) * ratio)
			i = getParent(tree,i)
		#~ else:
			#~ print getName(tree,l),"distance",distanceFrom(tree,root,l)



def constructSupportFromBootstrapTrees(tree,setoftrees):
  support = {}
  leaves = getLeavesNames(tree)
  pos = {}
  for i in range(len(leaves)):
    pos[leaves[i]] = i
  bipartitions = {}
  def complement(seq):
    result = []
    for s in seq.split("_")[0]:
      if s == "1":
	result.append("0")
      else:
	result.append("1")
    return "".join(result)
  def seq(leafset,node):
    result = ["0"]*len(leaves)
    for l in leafset:
      result[pos[l]] = "1"
    return "".join(result)
  def constructBipartAndReturnLeaves(tree,node):
    if isLeaf(tree,node):
      return [getName(tree,node)]
    else:
      c = getChildren(tree,node)
      tab0 = constructBipartAndReturnLeaves(tree,c[0])
      tab1 = constructBipartAndReturnLeaves(tree,c[1])
      result = tab0+tab1
      if len(c) > 2:
	tab2 = constructBipartAndReturnLeaves(tree,c[2])
	result = result + tab2
      if not isRoot(tree,node):
	support[node] = 0
	s = seq(tab0+tab1,node)
	bipartitions[s] = node
	bipartitions[complement(s)] = node
      return result

  root = getRoot(tree)
  constructBipartAndReturnLeaves(tree,root)

  def testBipartAndReturnLeaves(tree,node):
    #print "boot"
    if isLeaf(tree,node):
      return [getName(tree,node)]
    else:
      #print "nonleafboot"
      c = getChildren(tree,node)
      tab0 = testBipartAndReturnLeaves(tree,c[0])
      tab1 = testBipartAndReturnLeaves(tree,c[1])
      result = tab0+tab1
      if len(c) > 2:
	tab2 = testBipartAndReturnLeaves(tree,c[2])
	result = result + tab2
      if not isRoot(tree,node):
	s = seq(tab0+tab1,node)
	#print s
	if bipartitions.has_key(s):
	  #print "bip trouve"
	  support[bipartitions[s]] = support[bipartitions[s]] + 1
	#if bipartitions.has_key(complement(s)):
	  #support[bipartitions[complement(s)] = support[bipartitions[complement(s)]] + 1
      return result

  for t in setoftrees:
    root = getRoot(t)
    testBipartAndReturnLeaves(t,root)
  if len(setoftrees) > 0:
    for k in support.keys():
      writeBootstrap(tree,k,support[k]/float(len(setoftrees)))
  #root = getRoot(tree)
  #print writeTree(tree,root,False)
    

def RF(arbre1,arbre2):
  
  root1 = getRoot(arbre1)
  root2 = getRoot(arbre2)

  nodes1 = getNodes(arbre1)
  nodes2 = getNodes(arbre2)

  clades1 = []
  for n in nodes1:
	  leaves = getLeaves(arbre1,n)
	  if len(leaves) > 1:
		  clade = []
		  for l in leaves:
			  clade.append(getName(arbre1,l).split("|")[0].split("__")[0])
		  clade.sort()
		  clades1.append(clade)

  clades2 = []
  for n in nodes2:
	  leaves = getLeaves(arbre2,n)
	  if len(leaves) > 1:
		  clade = []
		  for l in leaves:
			  clade.append(getName(arbre2,l).split("|")[0].split("__")[0])
		  clade.sort()
		  clades2.append(clade)

  distance = 0
  for c in clades1:
	  if not c in clades2:
		  distance = distance + 1
		  #print 1,c
  for c in clades2:
	  if not c in clades1:
		  distance = distance + 1
		  #print 2,c

  return distance/2



def commonTriplets(arbre1,arbre2):
  result = 0
  triplets = {}
  
  root1 = getRoot(arbre1)
  leaves1 = getLeaves(arbre1,root1)
  for n1 in range(len(leaves1)):
    for n2 in range(n1+1,len(leaves1)):
      for n3 in range(n2+1,len(leaves1)):
	ids = [leaves1[n1],leaves1[n2],leaves1[n3]]
	ids.sort(lambda x,y: cmp(getName(arbre1,x),getName(arbre1,y)))
	names = [getName(arbre1,ids[0]),getName(arbre1,ids[1]),getName(arbre1,ids[2])]
	LCA12 = lastCommonAncestor(arbre1,ids[0],ids[1])
	LCA13 = lastCommonAncestor(arbre1,ids[0],ids[2])
	LCA23 = lastCommonAncestor(arbre1,ids[1],ids[2])
	#print LCA12,LCA13,LCA23
	if LCA12 == LCA13:
	  triplets['_'.join(names)] = 1
	if LCA12 == LCA23:
	  triplets['_'.join(names)] = 2
	if LCA13 == LCA23:
	  triplets['_'.join(names)] = 3
	#print names,triplets['_'.join(names)]
	  
  
  root2 = getRoot(arbre2)
  leaves2 = getLeaves(arbre2,root2)
  for n1 in range(len(leaves2)):
    for n2 in range(n1+1,len(leaves2)):
      for n3 in range(n2+1,len(leaves2)):
	#print n1,n2,n3,result
	ids = [leaves2[n1],leaves2[n2],leaves2[n3]]
	ids.sort(lambda x,y: cmp(getName(arbre2,x),getName(arbre2,y)))
	names = [getName(arbre2,ids[0]),getName(arbre2,ids[1]),getName(arbre2,ids[2])]
	if triplets.has_key('_'.join(names)):
	  LCA12 = lastCommonAncestor(arbre2,ids[0],ids[1])
	  LCA13 = lastCommonAncestor(arbre2,ids[0],ids[2])
	  LCA23 = lastCommonAncestor(arbre2,ids[1],ids[2])
	  if LCA12 == LCA13 and triplets['_'.join(names)] == 1:
	    #print names,"yes",triplets['_'.join(names)]
	    result = result + 1
	  elif LCA12 == LCA23 and triplets['_'.join(names)] == 2:
	    #print names,"yes",triplets['_'.join(names)]
	    result = result + 1
	  elif LCA13 == LCA23 and triplets['_'.join(names)] == 3:
	    #print names,"yes",triplets['_'.join(names)]
	    result = result + 1
	  #else:
	    #print names
	#else:
	  #print names,"not found"

  return result

		
# structure of the tree:
# 0: name, 1: parent, 2: tab of children, 3: length, 4: isdup, 5:species, 6:bootstrap 




#####################################################
#####################################################
#  Traversal of one tree 
# 
#####################################################
#####################################################


def readTree(treeseq):

    ###############################################
    ######### TREE READING ########################
    ###############################################
    tree = {"sequence":treeseq}
    id_node = 0
    nb_parenth = 0
    bppnumber = 0
    pile = []
    t = 0
    while t < len(treeseq):
        if treeseq[t] == "(":
            id_node = id_node + 1
	    nb_parenth = nb_parenth + 1
            tree[id_node]={}
            tree[id_node][0] = "N"+str(id_node)
            tree[id_node][1] = -1
            tree[id_node][2] = []
            tree[id_node][3] = 0
            tree[id_node][4] = ""
            tree[id_node][5] = ""
            tree[id_node][6] = ""
            tree[id_node][7] = -1
                        # [nom,pere,[enfants],longueur,annotD,dotannot,bootstrap,bppnumber]
#            print "ouverture",tree[id_node]
            if len(pile) > 0:
                tree[id_node][1] = pile[-1]
            pile.append(id_node)
            t = t + 1
        elif treeseq[t] == ")":
            t = t + 1
	    nb_parenth = nb_parenth - 1
	    tree[pile[-1]][7] = bppnumber
	    bppnumber = bppnumber + 1
	    #~ print nb_parenth,"(-1)",treeseq[t:t+80]
 
            if treeseq[t] == "@":
                t = t + 1
                tree["ancestor"] = pile[-1]
    
            while (treeseq[t] != ":" and
			treeseq[t] != ";" and
			treeseq[t] != "[" and
			treeseq[t] != ")" and
			treeseq[t] != ","):
				tree[pile[-1]][6] = tree[pile[-1]][6] + treeseq[t]
				t = t + 1

            if treeseq[t] == ":":
                debut = t + 1
                while treeseq[t] != "," and treeseq[t]!=")" and treeseq[t] != "[" and treeseq[t] != ";":
                    t = t + 1
                longueur = float(treeseq[debut:t])
                tree[pile[-1]][3] = longueur
                while treeseq[t] != "," and treeseq[t] != ")" and treeseq[t] != "[" and treeseq[t] != ";":
                    t = t + 1

                    
            if treeseq[t] == "[":
                debut = t + 1
                t = debut + treeseq[debut:].find("]")
                chaine = treeseq[debut:t]
                mots = chaine.split(":")
                for m in mots:
                    if m == "D=Y" or m == "D=T" or m == "Ev=GDup":
                        tree[pile[-1]][4] = "D"
                    if m[:2] == "S=":
                        tree[pile[-1]][5] = m[2:]
                    if m[:2] == "B=":
		        tree[pile[-1]][6] = m[2:]
		    if m[:3] == "ND=":
		        tree[pile[-1]]["ND"] = m[3:]
		    if isfloat(m):
		      tree[pile[-1]][6] = float(m)
                t = t + 1

            if treeseq[t] == ":":
                debut = t + 1
                while treeseq[t] != "," and treeseq[t]!=")" and treeseq[t] != "[" and treeseq[t] != ";":
                    t = t + 1
                longueur = float(treeseq[debut:t])
                tree[pile[-1]][3] = longueur
                while treeseq[t] != "," and treeseq[t] != ")" and treeseq[t] != "[" and treeseq[t] != ";":
                    t = t + 1


            del pile[-1]
            
            if treeseq[t] == ";":
                t = len(treeseq)
                
        elif treeseq[t] == ";":
            t = len(treeseq)
            
        elif treeseq[t]==",":
            t = t + 1
	    
        elif treeseq[t]==" ":
            t = t + 1
            
        else:  # nom d'une feuille
	    #print "nom_de_feuille"
            id_node = id_node + 1
            tree[id_node] = {}
            tree[id_node][1] = -1
            tree[id_node][2] = []
            tree[id_node][3] = 0
            tree[id_node][4] = ""
            tree[id_node][5] = ""
            tree[id_node][6] = ""
            tree[id_node][7] = bppnumber
	    bppnumber = bppnumber + 1
            if len(pile)>0:
                tree[id_node][1]=pile[-1]
            pile.append(id_node)
            debut = t
            while (treeseq[t]!="," and
                   treeseq[t]!=")" and
                   treeseq[t]!=":" and
                   treeseq[t]!=";" and  
                   treeseq[t]!="\n" and
                   treeseq[t] != "["):
                t=t+1
            nom = treeseq[debut:t].strip()
            tree[pile[-1]][0] = nom
	    #~ print nom
            
            if treeseq[t]==":":
                debut = t + 1
                while treeseq[t]!="," and treeseq[t]!=")" and treeseq[t] != "[" and treeseq[t] != ";":
                    t = t + 1
                longueur = float(treeseq[debut:t])
                tree[id_node][3] = longueur
               
            #print "fin nom"
            if treeseq[t] == "[":
                debut = t + 1
                t = debut + treeseq[debut:].find("]")
                chaine = treeseq[debut:t]
                #print chaine
                mots = chaine.split(":")
                for m in mots:
                    if m[:2] == "S=":
                        tree[pile[-1]][5] = m[2:]
                    if m[:3] == "ND=":
		        tree[pile[-1]]["ND"] = m[3:]    
                t = t + 1
                
            if treeseq[t]==":":
                debut = t + 1
                while treeseq[t]!="," and treeseq[t]!=")" and treeseq[t] != "[" and treeseq[t] != ";":
                    t = t + 1
                longueur = float(treeseq[debut:t])
                tree[id_node][3] = longueur
               
            del pile[-1]
    #print tree
    # remplissage des enfants
    nodes = list(getNodes(tree))
    for node in nodes:
        if not isRoot(tree,node):
	    pere = getParent(tree,node)
            addChild(tree,pere,node)
        
    return tree



