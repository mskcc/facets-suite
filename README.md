## Helper scripts that wrap the FACETS CRAN Library
![Facets is Great](https://i.imgur.com/vP2tUvM.jpg "Aww Yeah FACETS Time!")
```
usage: facets doFacets [-h] [-c CVAL] [-s SNP_NBHD] [-n NDEPTH] [-m MIN_NHET]
                       [-pc PURITY_CVAL] [-ps PURITY_SNP_NBHD]
                       [-pn PURITY_NDEPTH] [-pm PURITY_MIN_NHET] [-d DIPLOGR]
                       [-g GENOME] [-f COUNTS_FILE] [-t TAG] [-D DIRECTORY]
                       [-r R_LIB] [-C SINGLE_CHROM] [-G GGPLOT2] [--seed SEED]

optional arguments:
  -h, --help            show this help message and exit
  -c CVAL, --cval CVAL  critical value for segmentation
  -s SNP_NBHD, --snp_nbhd SNP_NBHD
                        window size
  -n NDEPTH, --ndepth NDEPTH
                        threshold for depth in the normal sample
  -m MIN_NHET, --min_nhet MIN_NHET
                        minimum number of heterozygote snps in a segment used
                        for bivariate t-statistic during clustering of
                        segments
  -pc PURITY_CVAL, --purity_cval PURITY_CVAL
                        critical value for segmentation
  -ps PURITY_SNP_NBHD, --purity_snp_nbhd PURITY_SNP_NBHD
                        window size
  -pn PURITY_NDEPTH, --purity_ndepth PURITY_NDEPTH
                        threshold for depth in the normal sample
  -pm PURITY_MIN_NHET, --purity_min_nhet PURITY_MIN_NHET
                        minimum number of heterozygote snps in a segment used
                        for bivariate t-statistic during clustering of
                        segments
  -d DIPLOGR, --dipLogR DIPLOGR
                        diploid log ratio
  -g GENOME, --genome GENOME
                        Genome of counts file
  -f COUNTS_FILE, --counts_file COUNTS_FILE
                        paired Counts File
  -t TAG, --TAG TAG     output prefix
  -D DIRECTORY, --directory DIRECTORY
                        output prefix
  -r R_LIB, --R_lib R_LIB
                        Which version of FACETs to load into R
  -C SINGLE_CHROM, --single_chrom SINGLE_CHROM
                        Perform analysis on single chromosome
  -G GGPLOT2, --ggplot2 GGPLOT2
                        Plots using ggplot2
  --seed SEED           Set the seed for reproducability


```

```
usage: facets mafAnno [-h] [-m MAF] [-f FACETS_FILES] [-o OUT_MAF]

optional arguments:
  -h, --help            show this help message and exit
  -m MAF, --maf MAF     file name of maf file to be annotated.
  -f FACETS_FILES, --facets_files FACETS_FILES
                        Mapping of "Tumor_Sample_Barcode" from maf and
                        "Rdata_filename" from FACETS (tab-delimited with
                        header)
  -o OUT_MAF, --out_maf OUT_MAF
                        file name of CN annotated maf.
```

```
usage: facets normDepth [-h] [-f FILE]

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Filename of counts file to be normalized.

```

```
usage: facets geneLevel [-h] [-f FILENAMES] [-o OUTFILE]

optional arguments:
  -h, --help            show this help message and exit
  -f FILENAMES, --filenames FILENAMES
                        list of filenames to be processed.
  -o OUTFILE, --outfile OUTFILE
                        Output filename.
```

```
usage: facets mergeTN [-h] [-t TUMOR_COUNTS] [-n NORMAL_COUNTS] [-o OUTFILE]

optional arguments:
  -h, --help            show this help message and exit
  -t TUMOR_COUNTS, --tumor_counts TUMOR_COUNTS
                        Tumor counts file name
  -n NORMAL_COUNTS, --normal_counts NORMAL_COUNTS
                        Normal counts file name
  -o OUTFILE, --outfile OUTFILE
                        output file (gzipped)
```
