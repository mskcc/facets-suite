# facetsSuite
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build status](https://travis-ci.org/mskcc/facets-suite.svg?branch=Rpackagev2)](https://travis-ci.org/mskcc/facets-suite)
[![Coverage status](https://codecov.io/gh/mskcc/facets-suite/branch/Rpackagev2/graph/badge.svg)](https://codecov.io/github/mskcc/facets-suite?branch=Rpackagev2)

facetsSuite is an R package with functions to run [FACETS](https://github.com/mskcc/facets)—an allele-specific copy-number caller for paired tumor-normal DNA-sequencing data from genome-wide and targeted assays. facetSuite both wraps the code to execute the FACETS algorithm itself as well as performs _post-hoc_ analyses on the resulting data. This package was developed by members of the [Taylor lab](https://www.mskcc.org/research-areas/labs/barry-taylor) and the Computational Sciences group within the [Center for Molecular Oncology at Memorial Sloan Kettering Cancer Center](https://www.mskcc.org/research-programs/molecular-oncology).

## Installation

You can install facetsSuite in R from this repository with:

``` r
devtools::install_github("mskcc/facets-suite", ref = "Rpackagev2")
```

Also follow the [instructions for installing FACETS](https://github.com/mskcc/facets).

_Note: For the wrapper script `snp-pileup-wrapper.R` you need to specift the variable `snp_pileup_path` in the script to point to the installation path of snp-pileup._

## Usage

### R functions

The R functions in this package are documented and their description and usage is available in R by doing:
```r
?facetsSuite::function_name
```

### Wrapper scripts

Most use of this package can be done from the command line using three wrapper scripts:
- `snp-pileup-wrapper.R`:\
    This wraps the `snp-pileup` C++ script that genotypes sites across the genome in both normal and tumor samples. The output from this is the input to FACETS. Most default input arguments are appropriate regardless of usage, but `--max-depth` may need adjustment depending on the overall depth of the samples used.

- `run-facets-wrapper.R`:\
    This wrapper takes above SNP "pileup" as input and executes the FACETS algorithm. The ouputs are in the form of Rdata objects, TXT files, and PNGs of the samples overall copy-number profile. The wrapper allows for running FACETS in a two-pass mode, where first a "purity" run estimates the overall segmentation profile, sample purity and ploidy, and subsequently the dipLogR value from this run seeds a "high-sensitivity" run which may detect more focal events. To run in the two-pass mode, specify the input arguments prefixed by `purity`.

- `annotate-maf-wrapper.R`:\
    This script can be run for one or more samples with  somatic mutation calls in MAF format to estimate their cancer-cell fractions (CCFs).

All three wrappers use [argparse](https://github.com/trevorld/r-argparse) for argument handling and can thus be run with `--help` to see the available input arguments.
