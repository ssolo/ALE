# MaxTiC: Fast ranking of a phylogenetic tree by Maximum Time Consistency with lateral gene transfers

----------------------------------------------------------------------------------------------------------------------------------------------

## GENERAL INFORMATION

Python code written by Eric Tannier, Inria.

Using ideas, comments, suggestions from Cédric Chauve, Akbar Rafiey, Adrian A. Davin, Celine Scornavacca, Philippe Veber, Bastien Boussau, Gergely J Szöllosi, Vincent Daubin.

Software distributed under the [cecill licence](http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt), rights and permissions are described here:
http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

Bug reports, suggestions or help requests for usage of MaxTic should be sent to eric.tannier@inria.fr

If you use this software please cite
[MaxTiC: Fast ranking of a phylogenetic tree by Maximum Time Consistency with lateral gene transfers, Biorxiv doi.org/10.1101/127548](http://www.biorxiv.org/content/early/2017/04/14/127548)

----------------------------------------------------------------------------------------------------------------------------------------------


## INSTALLATION

You should have python 2.7 installed and have downloaded the two following Python files:

[MaxTic.py](https://github.com/ssolo/ALE/blob/master/misc/MaxTiC.py)
[script_tree.py](https://github.com/ssolo/ALE/blob/master/misc/script_tree.py)

Those two files **must** be placed in the same directory.

## USAGE
From the directory containing the two Python files type:

```
python MaxTiC.py species_tree_file constraints_file [ls=LOCAL_SEARCH]
```

Where:
* *species_tree_file*
 is the name of a file containing a phylogenetic tree in Newick format, with all internal nodes labeled in the usual bootstrap field. For example an appropriate file would contain:
`(((((E874:71176.953359,E876:71176.953359)147:91861.1927533,(E882:133963.416166,E884:133963.416166)148:29074.7299467)173:200567.267366,(((E906:81440.5931371,(E909:8211.56159025,E910:8211.56159025)149:73229.0315469)174:56454.5197908,E914:137895.112928)187:177491.76707,(E925:297995.062188,(E938:185198.008495,E944:185198.008495)150:112797.053694)175:17391.8178093)188:48218.5334805)197:1490.09088669,(E955:197786.332581,E957:197786.332581)151:167309.171784)203:445.96525336,(E967:274500.929682,((E982:139138.6704,E987:139138.6704)153:83916.622881,E990:223055.293281)154:51445.636401)177:91040.539936)208;`

Ultrametricity is not required, but will be used if provided, in order to compare the output and input rankings in terms of Kendall distance.

* *constraints_file*
 is the name of a file which contains a list of constraints between internal nodes of the species tree. One line of this file is a constraint, like in this example:
```
gene_family_1,147,149,0.09
gene_family_2,197,187,0.12
gene_family_2,187,188,0.9
```
 * the first field is a gene family identifier, which is not used by the program but has been useful for our results analyses
 * the second and third fields are labels of internal nodes of the species tree, telling that in a ranking, 147 should be older than 149, 197 should be older than 187, 187 should be older than 188
 * the last field is a weight associated with the constraint, which is supposed to be taken as a confidence score you can put on this constraint.

We construct the constraints file from the output of the software ALE, and the script to do so is available on demand. For each transfer detected by ALE, we report that the father of the donor branch should be older than the child of the receptor branch. However constraints can be constructed from any software detecting transfers or any type of data yielding relative time constraints between nodes.

* *ls=LOCAL_SEARCH*
 is an optional parameter that can improve your solution. LOCAL_SEARCH is the time in seconds during which you want to try improvements by a local search. Typically we made most analyses with ls=180.

---------------------------------------------------------------------------------------------------------------------------------------------

## OUTPUT

There are 4 kinds of output, the standard screen output and 3 additional files. The typical standard output is as follows:
```
tree with  108 internal nodes
8378.96 total weight of constraints from transfers
726.31 uninformative ( 726.31 to a descendant, 0 to a leaf 0 to an ancestor 0 to itself 7%)
794.43 trivially conflicting constraints (9%) (descendant to ancestor or trivial cycle)
value of the order given by the input tree 1037.14 (12.3779084755%)
value of the greedy heuristic 1280.95 (15.2876968025%)
value of the mixing heuristic: 974.38 (11.6288895042%)
best order is the mixing heuristic
((E23:20,(E25:6,E26:6)112:14)157:88,(((E80:89,((E94:65,((E96:17,E99:17)152:27,E110:44)176:21)189:20,(E124:67,E131:67)109:18)198:4)204:12,((E144:90,((E164:29,(E166:28,E169:28)110:1)155:10,E174:39)178:51)190:10,((((E211:37,(E220:33,(E231:26,E235:26)111:7)156:4)179:38,(E248:34,(E255:32,E260:32)113:2)158:41)180:11,(E285:74,E291:74)114:12)191:8,(((E299:30,E304:30)115:6,(E314:31,E320:31)116:5)159:22,(E326:22,E327:22)117:36)160:36)199:6)205:1)209:6,(((((E354:61,E356:61)118:26,((E367:25,E372:25)119:57,E383:82)120:5)161:17,(((E390:46,(E391:1,E392:1)121:45)162:46,(((E409:52,E412:52)122:12,(E436:51,E439:51)123:13)163:7,(E460:13,E461:13)124:58)164:21)181:7,E470:99)192:5)200:1,((((((E479:2,E480:2)125:48,(((E499:38,E501:38)126:10,E513:48)127:1,E515:49)128:1)165:22,E525:72)182:1,(((E533:27,E544:27)129:20,E555:47)130:13,E561:60)131:13)193:8,E574:81)201:21,(E586:63,E600:63)132:39)206:3)210:1,(((((E618:9,E624:9)133:9,E641:18)134:3,E647:21)135:33,(E655:41,(E658:16,E662:16)136:25)166:13)183:49,((((((E668:57,E674:57)137:13,(E677:5,E678:5)138:65)167:9,((E693:53,E703:53)139:24,((E709:4,E711:4)140:58,(E722:7,E725:7)141:55)168:15)184:2)194:12,(E735:88,(E745:56,(E753:11,E757:11)142:45)169:32)185:3)195:2,(((E772:69,((E779:15,(E780:14,E781:14)143:1)170:4,E787:19)186:50)196:7,((E798:40,(E800:35,E802:35)144:5)171:28,(E812:12,E813:12)145:56)172:8)202:8,(E831:8,E832:8)146:76)207:9)211:5,(((((E874:10,E876:10)147:35,(E882:43,E884:43)148:2)173:50,(((E906:23,(E909:3,E910:3)149:20)174:1,E914:24)187:59,(E925:78,(E938:55,E944:55)150:23)175:5)188:12)197:1,(E955:59,E957:59)151:37)203:1,(E967:80,((E982:42,E987:42)153:24,E990:66)154:14)177:17)208:1)212:5)213:3)214:1)215:1)216;
Similarity of the best order compared with the input order 0.94713918036
attempting a local search from the best found order, please wait 180.0 seconds
better solution 973.86
better solution 973.61
better solution 973.28
better solution 973.19
better solution 972.38
better solution 972.29
better solution 972.2
better solution 971.39
better solution 970.55
better solution 970.51
better solution 970.38
better solution 970.34
after local search 970.34 rejected
best found solution 970.34 (11.5806734965%)
((E23:16,(E25:7,E26:7)112:9)157:92,(((E80:89,((E94:65,((E96:15,E99:15)152:31,E110:46)176:19)189:20,(E124:67,E131:67)109:18)198:4)204:12,((E144:91,((E164:30,(E166:29,E169:29)110:1)155:9,E174:39)178:52)190:9,((((E211:37,(E220:34,(E231:28,E235:28)111:6)156:3)179:38,(E248:31,(E255:24,E260:24)113:7)158:44)180:11,(E285:74,E291:74)114:12)191:8,(((E299:32,E304:32)115:4,(E314:33,E320:33)116:3)159:22,(E326:23,E327:23)117:35)160:36)199:6)205:1)209:6,(((((E354:61,E356:61)118:26,((E367:27,E372:27)119:55,E383:82)120:5)161:17,(((E390:48,(E391:1,E392:1)121:47)162:44,(((E409:52,E412:52)122:12,(E436:51,E439:51)123:13)163:7,(E460:13,E461:13)124:58)164:21)181:7,E470:99)192:5)200:1,((((((E479:2,E480:2)125:48,(((E499:38,E501:38)126:3,E513:41)127:4,E515:45)128:5)165:22,E525:72)182:1,(((E533:18,E544:18)129:31,E555:49)130:11,E561:60)131:13)193:8,E574:81)201:21,(E586:63,E600:63)132:39)206:3)210:1,(((((E618:10,E624:10)133:10,E641:20)134:2,E647:22)135:32,(E655:42,(E658:17,E662:17)136:25)166:12)183:49,((((((E668:57,E674:57)137:13,(E677:6,E678:6)138:64)167:9,((E693:53,E703:53)139:24,((E709:4,E711:4)140:58,(E722:8,E725:8)141:54)168:15)184:2)194:11,(E735:88,(E745:56,(E753:11,E757:11)142:45)169:32)185:2)195:3,(((E772:69,((E779:19,(E780:14,E781:14)143:5)170:2,E787:21)186:48)196:7,((E798:40,(E800:35,E802:35)144:5)171:28,(E812:12,E813:12)145:56)172:8)202:8,(E831:9,E832:9)146:75)207:9)211:5,(((((E874:5,E876:5)147:42,(E882:44,E884:44)148:3)173:48,(((E906:25,(E909:3,E910:3)149:22)174:1,E914:26)187:57,(E925:78,(E938:55,E944:55)150:23)175:5)188:12)197:1,(E955:59,E957:59)151:37)203:1,(E967:80,((E982:43,E987:43)153:23,E990:66)154:14)177:17)208:1)212:5)213:3)214:1)215:1)216;
['216', '215', '214', '210', '200', '213', '206', '209', '205', '192', '212', '208', '203', '197', '199', '211', '181', '190', '195', '204', '185', '161', '191', '198', '207', '188', '120', '201', '177', '194', '175', '184', '202', '180', '114', '193', '182', '164', '167', '196', '172', '109', '154', '189', '163', '132', '168', '118', '131', '151', '160', '137', '169', '150', '183', '139', '122', '123', '165', '130', '162', '173', '176', '128', '148', '153', '166', '127', '171', '178', '126', '179', '159', '144', '156', '116', '115', '158', '155', '110', '111', '119', '187', '174', '113', '117', '135', '186', '134', '170', '129', '136', '157', '152', '143', '124', '145', '142', '133', '146', '141', '112', '138', '147', '140', '149', '125', '121']
0.0210997459716  constraints in agreement with the best tree in conflict with input tree
Similarity of the order compared with the input order 0.947931102752
```



The total weight is the sum of the weights of all transfers you provided as input. The program counts and discards the uninformative ones (constraint from an ancestor to a descendant for example). Then a ranking of the species tree maximizing the cumulative weight of non conflicting constraints is sought. The value of the input tree is given, meaningful only if the input tree is ultrametric. Then two heuristics are run and the result of the best one is output (numbers and percentages are the cumulative weight of constraints contradicting the solution, here 11.62% is the best so far).
Similarity is the Kendall normalized similarity, that is, ( MAX_DIST - DIST ) / MAX_DIST, where DIST is the standard Kendall distance between orders of internal nodes, and MAX_DIST is the maximum possible such distance given the input species tree. So Kendall similarity is a number between 0 and 1, 1 meaning exactly the same ranked tree.
At the end the solution of the local search from the best heuristic solution is output, the ranked tree, the internal node order and the similarity with input tree.


Then the three output files are:
* `constraints_file_MT_output_filtered_list_of_weighted_informative_constraints`: input file expurged from noninformative constraints
* `constraints_file_MT_output_list_of_constraints_conflicting_with_best_order`: list of constraints conflicting with the solution (ranked species tree)
* `constraints_file_MT_output_partial_order`: list of constraints agreeing with the solution (ranked species tree)
