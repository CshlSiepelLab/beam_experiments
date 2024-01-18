#A bayesian approach to inferring anatomical tissue locations for internal nodes in a cancer phylogeny with labeled leaves.

To simulate data from an activated simulate conda environment:
```
./scripts/sim_wrapper.sh --design RANDOM --out test_sim --sites 10 --mutrate 0.5 --samples 10 --migration migration_matrix/true_migration_prob_matrix.csv
```

The output should provide all that is necessary as input to the inference method, which is the true tree and the tissue labels of the leaves. The barcode data is also provided in the form of a mutation matrix.
