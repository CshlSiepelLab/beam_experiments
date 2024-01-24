# A bayesian approach to inferring anatomical tissue locations for internal nodes in a cancer phylogeny with labeled leaves.

### Simulating data

To simulate data from an activated simulate conda environment:
```
./scripts/sim_wrapper.sh --design RANDOM --out test_sim --sites 10 --mutrate 1.0 --samples 50 --migration inputs/test_migration_prob_matrix.csv
```

The output should provide all that is necessary as input to the inference method, which is the true tree and the tissue labels of the leaves. The barcode data is also provided in the form of a mutation matrix.

### Formatting simulated data to input to BEAST2 FixedTreeAnalysis

Use simulate environment and input tree newick file and tsv of node to tissue mapping from simulated data: NOT FUNCTIONAL YET
```
python ./scripts/format_fixed_tree_from_sim.py examples/simulated_data/sim_results_test_sim/test_sim_true.nwk examples/simulated_data/sim_results_test_sim/test_sim_tissues.tsv
```

### Formatting simulated data to input to TideTree and running TideTree
```
python scripts/format_tidetree_sequences_sim_matrix.py examples/simulated_data/sim_results_ten_samples/ten_samples_indel_character_matrix.tsv
scripts/format_tidetree_xml_from_sim.sh
scripts/run_tidetree.sh --xml tidetree.xml
```

