# Comparing models with reseeding allowed vs not allowed for MMUS1467 CP01 in Serio et al. data

Without reseeding:
```
java -Xmx10g -jar /grid/siepel/home_norepl/staklins/beam/beam.jar -threads 5 -overwrite -working /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/pipelines/model_selection/beam_ns_template_no_reseeding_mmus1467_cp01.xml
```

With reseeding:
```
java -Xmx10g -jar /grid/siepel/home_norepl/staklins/beam/beam.jar -threads 5 -overwrite -working /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/pipelines/model_selection/beam_ns_template_one_rate_reseeding_mmus1467_cp01.xml
```

Then, compute the Bayes Factor from the marginal likelihoods estimated for each model.

