#!/usr/bin/env python3

import sys
from ete3 import Tree

newick = sys.argv[1]
id = sys.argv[2]
outdir = sys.argv[3]

# newick = '((((1[&location="P"]:27.097507928285754,4[&location="P"]:27.097507928285754)[&location="P"]:0.20708955492136383,(((((7[&location="P"]:23.091612927259938,(10[&location="M1"]:12.807363959492575,((17[&location="M1"]:10.442452997696131,(((27[&location="M4"]:2.723126636364711,(48[&location="M4"]:2.5419538299728837,49[&location="M4"]:2.5419538299728837)[&location="M4"]:0.1811728063918272)[&location="M4"]:0.11307166605351027,47[&location="M4"]:2.836198302418221)[&location="M4"]:0.37005061280697404,50[&location="M4"]:3.206248915225195)[&location="M4"]:7.236204082470936)[&location="M1"]:0.9600906706491052,(18[&location="M1"]:9.535251509159636,19[&location="M1"]:9.535251509159636)[&location="M1"]:1.8672921591856007)[&location="M1"]:1.4048202911473382)[&location="M1"]:10.284248967767363)[&location="P"]:0.704566121236887,11[&location="P"]:23.796179048496825)[&location="P"]:0.8709118726188159,((28[&location="M8"]:10.095578025062663,(29[&location="M8"]:6.408952206871983,30[&location="M8"]:6.408952206871983)[&location="M8"]:3.68662581819068)[&location="M8"]:0.21340281857226273,(39[&location="M8"]:8.505859473212853,40[&location="M8"]:8.505859473212853)[&location="M8"]:1.8031213704220725)[&location="M8"]:14.358110077480715)[&location="P"]:0.49619051254560986,8[&location="P"]:25.16328143366125)[&location="P"]:1.5048938875498372,((((12[&location="M3"]:18.613764328743173,(((20[&location="M3"]:15.919269258519288,31[&location="M3"]:15.919269258519288)[&location="M3"]:1.666026937751914,33[&location="M3"]:17.5852961962712)[&location="M3"]:0.09582823357022718,42[&location="M3"]:17.68112442984143)[&location="M3"]:0.9326398989017441)[&location="M3"]:0.9701634278179618,41[&location="M3"]:19.583927756561135)[&location="M3"]:0.9997127438963069,21[&location="M3"]:20.58364050045744)[&location="M3"]:1.7885260989039828,32[&location="M3"]:22.372166599361424)[&location="M3"]:4.296008721849663)[&location="P"]:0.6364221619960304)[&location="P"]:0.20288081084065368,3[&location="P"]:27.507478294047772)[&location="P"]:3.4802954829650012,(((2[&location="P"]:25.362466717194977,6[&location="P"]:25.362466717194977)[&location="P"]:0.26790503091293516,(9[&location="M9"]:3.8826053009235295,(13[&location="M9"]:2.291932102452617,14[&location="M9"]:2.291932102452617)[&location="M9"]:1.5906731984709124)[&location="M9"]:21.747766447184382)[&location="P"]:0.7410918496964456,(5[&location="P"]:25.833565864673986,(15[&location="M7"]:16.10741091174917,(((16[&location="M7"]:13.371539740825696,(22[&location="M7"]:10.914799234764939,((23[&location="M7"]:10.035414309594168,24[&location="M7"]:10.035414309594168)[&location="M7"]:0.6457552113611094,34[&location="M7"]:10.681169520955278)[&location="M7"]:0.2336297138096608)[&location="M7"]:2.4567405060607577)[&location="M7"]:0.7379092597908841,(((((((25[&location="M6"]:1.5156184062841176,26[&location="M6"]:1.5156184062841176)[&location="M6"]:1.0666857628850053,45[&location="M6"]:2.582304169169123)[&location="M6"]:1.0699327025335235,43[&location="M6"]:3.6522368717026463)[&location="M6"]:0.07693274196284206,44[&location="M6"]:3.7291696136654884)[&location="M6"]:0.09070922573797224,(36[&location="M6"]:3.080336187365207,46[&location="M6"]:3.080336187365207)[&location="M6"]:0.7395426520382538)[&location="M6"]:1.6695397728926658,37[&location="M6"]:5.489418612296126)[&location="M6"]:1.0446838467845323,38[&location="M6"]:6.534102459080659)[&location="M6"]:7.575346541535922)[&location="M7"]:1.4434675201193397,35[&location="M7"]:15.55291652073592)[&location="M7"]:0.5544943910132503)[&location="M7"]:9.726154952924816)[&location="P"]:0.5378977331303716)[&location="P"]:4.616310179208416)[&location="P"]:0.0;'
# id = 10000
# outdir = '.'

# Read the BEAST tree with tissue annotations
tree = Tree(newick, format=1)

# convert the annotation to name labels
tree_beast_annotation = tree.copy()
i = 0
for node in tree_beast_annotation.traverse():
    name = node.name
    if name.startswith("[&location="):
        name = f"node{i}{name}"
        i += 1
    name = name.replace('[&location="', "_")
    name = name.replace('"]', "")
    node.name = name

# remove annotations
tree_no_annotations = tree_beast_annotation.copy()
for node in tree_no_annotations.traverse():
    node.name = node.name.split("_")[0]

# get a tsv file of the annotations for the tips
with open(f"{outdir}/{id}_tip_tissues.tsv", "w") as f:
    for node in tree_beast_annotation.iter_leaves():
        name, tissue = node.name.split("_")
        f.write(f"{name}\t{tissue}\n")


# Write the trees to files
tree_no_annotations.write(outfile=f"{outdir}/{id}.newick", format=8)
tree_beast_annotation.write(outfile=f"{outdir}/{id}_beast.newick", format=8)
