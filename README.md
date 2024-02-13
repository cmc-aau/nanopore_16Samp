# Snakemake workflow: `cmc-aau/nanopore_16Samp`

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.18.2-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/cmc-aau/nanopore_16Samp/workflows/Tests/badge.svg)](https://github.com/cmc-aau/nanopore_16Samp/actions?query=branch%3Amain+workflow%3ATests)

Workflow to map any amplicon reads against a taxonomic database and produce an abundance table compatible with the [ampvis2](https://github.com/kasperskytte/ampvis2) R package.

## Usage
The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=cmc-aau/nanopore_16Samp).

First install `snakemake` and `snakedeploy`, then deploy the workflow in your project using fx

```
snakedeploy deploy-workflow https://github.com/cmc-aau/nanopore_16Samp . --tag v123
```

Replace `v123` with one of the versions listed under releases.

Then adjust the configuration file `config.yml` and run the workflow with `snakemake` according to your setup.

## References
https://doi.org/10.1016/j.watres.2023.119919
