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

_Note: For the wrapper script `snp-pileup-wrapper.R` you need to specify the variable `snp_pileup_path` in the script to point to the installation path of snp-pileup._

## Usage

### R functions

The R functions in this package are documented and their description and usage is available in R by doing:
```r
?facetsSuite::function_name
```

Central to most functionality in the package is the output from the `run_facets`, which runs the FACETS algorithm based on provided tumor-normal SNP pileup (i.e. genotyping). The output is a list object with the following named objects:
- `snps`: SNPs used for copy-number segmentation, where `het==1` indicates heterozygous loci.
- `segs`: Inferred copy-number segmentation.
– `purity`: Inferred sample purity, i.e. fraction of tumor cells of the total cellular population.
- `ploidy`: Inferred sample ploidy.
- `diplogr`: Inferred dipLogR, the sample-specific baseline corresponding to the diploid state.
- `alballogr`: Alternative dipLogR value(s) at which a balanced solution was found.
- `flags`: Warning flags from the naïve segmentation algorithm.
- `em_flags`: Warning flags from the expectation-maximization segmentation algorithm.
- `loglik`: Log-likelihood value of the fitted model.

Note that FACETS performs segmentation with two algorithms, the "naïve" base method and an expectation-maximization algorithm. The latter (columns suffixed `.em`) is used as a default for most of the functions in this package.

### Wrapper scripts

Most use of this package can be done from the command line using three wrapper scripts:
- `snp-pileup-wrapper.R`:\
    This wraps the `snp-pileup` C++ script that genotypes sites across the genome in both normal and tumor samples. The output from this is the input to FACETS. Most default input arguments are appropriate regardless of usage, but `--max-depth` may need adjustment depending on the overall depth of the samples used.\
    Example command:
    ```shell
    snp-pileup-wrapper.R \
        --snp-pileup-path <path to snp-pileup executable> \
        --vcf-file <path to SNP VCF> \
        --normal-bam normal.bam \
        --tumor-bam tumor.bam
    ```
    The input VCF file should contain polymorphic SNPs, so that FACETS can infer changes in allelic configuration at genomic loci from changes in allele ratios. [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) is a good source for this. By default, `snp-pileup` also estimates the read depth in the input BAM files every 50th base.

- `run-facets-wrapper.R`:\
    This wrapper takes above SNP "pileup" as input and executes the FACETS algorithm. The ouputs are in the form of Rdata objects, TXT files, and PNGs of the samples overall copy-number profile. The wrapper allows for running FACETS in a two-pass mode, where first a "purity" run estimates the overall segmentation profile, sample purity and ploidy, and subsequently the dipLogR value from this run seeds a "high-sensitivity" run which may detect more focal events. To run in the two-pass mode, specify the input arguments prefixed by `purity`. The cval (`--purity-cval` and `--cval`) parameters tune the segmentation coarseness.\
    Example command:
    ```shell
    run-facets-wrapper.R \
        --counts-file tumor_normal.snp_pileup.gz \
        --sample-id tumorID__normalID \
        --purity-cval 1000 --cval 500 \
        --everything
    ```
    The above command runs FACETS in the two-pass mode, first at cval 1000, then at cval 500 based on the sample-specific baseline found at the higher cval. The full suite of analysis and QC is run with the `--everything` flag. If no output directory is specified, a directory named `sample-id` is created.

- `annotate-maf-wrapper.R`:\
    This script estimates the cancer-cell fractions (CCFs) of somatic mutations using purity and ploidy estimates from FACETS. It requires a input MAF file and a mapping of sample names in the MAF file (column `Tumor_Sample_Barcode`) to FACETS output RDS files (i.e. file paths). Alternatively, it can be run in a single-sample mode by pointing direct to the RDS and providing a MAF file with only mutation calls for the given sample.\
    Example command:
    ```shell
    annotate-maf-wrapper.R \
        --maf-file somatic_mutations.maf
        --facets-output <path to facets_output.rds>
    ```
    Or run with a mapping file as input (`--sample-mapping`), in the following format:
    ```shell
    > cat sample_map.txt
    sample      file
    SampleA     SampleA_facets.rds
    SampleB     SampleB_facets.rds
    ...         ...
    ```


All three wrappers use [argparse](https://github.com/trevorld/r-argparse) for argument handling and can thus be run with `--help` to see the all input arguments.
