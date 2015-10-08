## Helper scripts that wrap the FACETS CRAN Library
![Facets is Great](https://i.imgur.com/vP2tUvM.jpg "Aww Yeah FACETS Time!")
```
usage: facets [-h] {call,mergeTN,run,norm,maf} ...

run FACETS analysis

positional arguments:
  {call,mergeTN,run,norm,maf}
                        sub-command help
    run                 run FACETS from merged counts files
    mergeTN             merge tumor and normal outputs
    norm                run FACETS from merged counts files
    maf                 use FACETS results to annotate a maf file
    call                extract gene copy number status from FACETS results

optional arguments:
  -h, --help            show this help message and exit
  ```
  
```
usage: facets maf [-h] maf_file facets_sample_file output_file

positional arguments:
  maf_file            maf file for anntation, Tumor_Sample_Barcodes must
                      correspond to FACETS sample file
  facets_sample_file  file mapping each "Tumor_Sample_Barcode" from the maf
                      file to an Rdata_filename FACETS result (tab-delimited
                      with header)
  output_file         annotated maf file output filename
```

```
usage: facets run [-h] [--lib-version LIB_VERSION]
                  count_file tag output_directory

positional arguments:
  count_file            merged tumor/normal count file
  tag                   merged tumor/normal count file
  output_directory      merged tumor/normal count file

optional arguments:
  -h, --help            show this help message and exit
  --lib-version LIB_VERSION
                        path to facets library, if you don't want to use
                        default
```

```
usage: facets call [-h] [-o OUTPUT_FILE]
                   facets_cncf_files [facets_cncf_files ...]

positional arguments:
  facets_cncf_files     FACETS cncf.txt files
```
