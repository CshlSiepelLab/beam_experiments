#!/use/bin/env python3

from ete3 import Tree

machina_nwk = sys.argv[1]
# machina_nwk="data/longer_10million_mcmc_unsymmetrical_machina_m5_sims_compare_beast_machina_fixedtreeanalysis_default_2_12_24/machina_m5_sim_data/seed49/machina_tree_all_tissue_labels.nwk"

outputfile = machina_nwk.split(".")[0] + ".tree"

machina_tree = Tree(machina_nwk, format=8)

for node in machina_tree.traverse():
    if "_" in node.name:
        name_parts = node.name.split("_")
        metadata = name_parts[-1]
        node.add_feature("location", metadata)
        node.name = name_parts[0]

nex = machina_tree.write(features=["location"], format=9)

nex_new = nex.replace("&&NHX:", "&")

with open(outputfile, 'w') as file:
    file.write(f"#NEXUS\nBegin trees;\ntree TREE1={nex_new}\nEnd;")
