# A bayesian approach to inferring anatomical tissue locations for internal nodes in a cancer phylogeny with labeled leaves.

### Simulating data

To simulate data from an activated simulate conda environment:
```
./scripts/sim_wrapper.sh --design RANDOM --out test_sim --sites 10 --mutrate 1.0 --samples 50 --migrationrate 1.0 --migration inputs/test_migration_prob_matrix.csv
```

The output should provide all that is necessary as input to the inference method, which is the true tree and the tissue labels of the leaves. The barcode data is also provided in the form of a mutation matrix.

### Formatting simulated data to input to BEAST2 FixedTreeAnalysis

Use simulate environment and input tree newick file and tsv of node to tissue mapping from simulated data: NOT FUNCTIONAL YET
```
python ./scripts/format_fixed_tree_from_sim.py examples/simulated_data/sim_results_test_sim/test_sim_true.nwk examples/simulated_data/sim_results_test_sim/test_sim_tissues.tsv
```

### Formatting simulated data to input to TideTree and running TideTree

Takes in a mutation matrix and formats xml sequences section:
```
python scripts/format_tidetree_sequences_sim_matrix.py examples/simulated_data/sim_results_ten_samples/ten_samples_indel_character_matrix.tsv
```

Takes in xml sequences formatted and experimental parameters to fill in full xml template for TideTree:
```
scripts/format_tidetree_xml_from_sim.sh --seqs inputs/tidetree_seqs_example.xml --total 54 --edit 36 --chain 1000000
```

Runs TideTree from a fully formatted xml:
```
scripts/run_tidetree.sh --xml tidetree.xml
```

### Running FixedTreeAnalysis (or SetTreeAnalysisfrom TideTree output)

This can be done by following those respective tutorials which specify to import `.tree` file and `.dat` file for tissue locations into BEAUTi where substitution model parameters and MCMC options can be specified through the GUI.

This can also be done by directly editing the xml file as done in the pipelines in this repo.

### Resolving branch length estimated for MACHINA simulated data following PathFinder method

`python ./scripts/pathfinder/pathfinder.py input.fas --primary P -o Example_output`

### Running Bayes factor comparison for two XML files for FixedTreeAnalysis

```
scripts/bayes_factor_nested_sampling_from_xmls.sh --xml1 results/no_bsvss_compare_beast_machina_fixedtreeanalysis_3_4_24/no_bsvss_symmmetricalfalse_machina_m5_sims_compare_beast_machina_fixedtreeanalysis_3_4_24/machina_m5_sim_data/seed3/T_seed3_unlabeled_true_tree_final_input_xml.xml --xml2 results/no_bsvss_compare_beast_machina_fixedtreeanalysis_3_4_24/no_bsvss_symmmetricaltrue_machina_m5_sims_compare_beast_machina_fixedtreeanalysis_3_4_24/machina_m5_sim_data/seed3/T_seed3_unlabeled_true_tree_final_input_xml.xml --dir compare_asym_sym_bayes_factor_m5_seed3_3_5_24
```
or the more general form:
```
bayes_factor_nested_sampling_from_xmls.sh --xml1 <xml filepath (str)> --xml2 <xml filepath (str)> --dir <working directory path (str)>
```
