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

Need to provide xml file for a normal MCMC run of two models and then this script will reformat the xml file to run the Nested Sampling BEAST2 package to obtain marginal likelihood values with standard deviations that can be use dto get Bayes factor and compared to random chance differences:

```
scripts/bayes_factor_nested_sampling_from_xmls.sh --xml1 results/no_bsvss_compare_beast_machina_fixedtreeanalysis_3_4_24/no_bsvss_symmmetricalfalse_machina_m5_sims_compare_beast_machina_fixedtreeanalysis_3_4_24/machina_m5_sim_data/seed3/T_seed3_unlabeled_true_tree_final_input_xml.xml --xml2 results/no_bsvss_compare_beast_machina_fixedtreeanalysis_3_4_24/no_bsvss_symmmetricaltrue_machina_m5_sims_compare_beast_machina_fixedtreeanalysis_3_4_24/machina_m5_sim_data/seed3/T_seed3_unlabeled_true_tree_final_input_xml.xml --dir compare_asym_sym_bayes_factor_m5_seed3_3_5_24
```
or the more general form:
```
bayes_factor_nested_sampling_from_xmls.sh --xml1 <xml filepath (str)> --xml2 <xml filepath (str)> --dir <working directory path (str)>
```

It is also possible to seperate the steps for calculating the marginal likelihood of each model and then computing Bayes factor for both results as follows:
```
scripts/nested_sampling_marginal_likelihood_from_xml.sh --xml results/no_bsvss_compare_beast_machina_fixedtreeanalysis_3_4_24/no_bsvss_symmmetricalfalse_machina_m5_sims_compare_beast_machina_fixedtreeanalysis_3_4_24/machina_m5_sim_data/seed3/T_seed3_unlabeled_true_tree_final_input_xml.xml  --dir compare_asym_sym_bayes_factor_m5_seed3_3_5_24/model1

scripts/nested_sampling_marginal_likelihood_from_xml.sh --xml results/no_bsvss_compare_beast_machina_fixedtreeanalysis_3_4_24/no_bsvss_symmmetricaltrue_machina_m5_sims_compare_beast_machina_fixedtreeanalysis_3_4_24/machina_m5_sim_data/seed3/T_seed3_unlabeled_true_tree_final_input_xml.xml --dir compare_asym_sym_bayes_factor_m5_seed3_3_5_24/model2

scripts/bayes_factor_from_marginals.sh --ml1 compare_asym_sym_bayes_factor_m5_seed3_3_5_24/model1/xml1/xml1_marginal_likelihood_run.txt --ml2 /home/staklins/bayesian_phylogenetic_metastasis/compare_asym_sym_bayes_factor_m5_seed3_3_5_24/model2/xml1/xml1_marginal_likelihood_run.txt --dir compare_asym_sym_bayes_factor_m5_seed3_3_5_24/bayes_factor
```

Any run of Nesdted Sampling is dependent on 2 parameters, the number of active particles and the subchain length which need to be tuned appropriately according to Nested Sampling wiki. The general premise is that more active particles decreases the standard deviation of each marginal likelihood estimate, which then lowers the threshold for a Bayes factor difference to not just be caused by random chance. The sub chain length is less well defined but has to do with independent sampling of the next points from the current position and needs to be large enough to where further increases in length no longer make a difference in the estimate. This can be roughly estimated in the script `scripts/determine_subchainlength_nested_sampling.sh`.

### Running new BEAST2.7 metastabayes package models for reduced paramaterization and/or joint inference with TideTree

Need to obtain the metastabayes.jar from the metastabayes repo, which if placed up one directory from this repo, then can be run as follows with the option of adding the `-working` flag before the xml file name to produce output in the same directory as the xml file and the `-overwrite` flag to force BEAST to run a new analysis and overwrite previous results from the same xml in the same directory:
```
java -jar ../metastabayes/metastabayes.jar ../metastabayes/examples/oneRate_machina_m5_seed3.xml
```

