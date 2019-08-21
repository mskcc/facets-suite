#!/usr/bin/env Rscript

# Generate gene-level TCN/LCN and CNA call, given cncf.txt files from a doFacets run

## Require *hisens.cncf.txt

library(argparse)
library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(jsonlite)

"%ni%" = Negate("%in%")

write.text <- function (...) {
    write.table(..., quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}

getSDIR <- function(){
    args=commandArgs(trailing=FALSE)
    TAG="--file="
    path_idx=grep(TAG,args)
    SDIR=dirname(substr(args[path_idx],nchar(TAG)+1,nchar(args[path_idx])))
    if(length(SDIR)==0) {
        return(getwd())
    } else {
        return(SDIR)
    }
}

# Get IMPACT341 loci and gene names
IMPACT341_targets = function() {
    IMPACT341_targets = suppressWarnings(fread(paste0('grep -v "^@" ',getSDIR(),"/data/impact341_targets.ilist")))
    setnames(IMPACT341_targets, c("chr", "start", "end", "strand", "name"))
    setkey(IMPACT341_targets, chr, start, end)
    IMPACT341_targets
} 

# Get IMPACT410 loci and gene names
IMPACT410_targets = function() {
    IMPACT410_targets <- suppressWarnings(fread(paste0('grep -v "^@" ',getSDIR(),"/data/impact410_targets.ilist")))
    setnames(IMPACT410_targets, c("chr", "start", "end", "strand", "name"))
    setkey(IMPACT410_targets, chr, start, end)
    IMPACT410_targets
}

# Get IMPACT468 loci and gene names
IMPACT468_targets = function() {
    IMPACT468_targets = suppressWarnings(fread(paste0('grep -v "^@" ',getSDIR(),"/data/impact468_targets.ilist")))
    setnames(IMPACT468_targets, c("chr", "start", "end", "strand", "name"))
    setkey(IMPACT468_targets, chr, start, end)
    IMPACT468_targets
}

# Get exome-wide targets
exome_targets = function() {
    exome_targets = suppressWarnings(fread(paste0('grep -v "^@" ',getSDIR(),"/data/Homo_sapiens.GRCh37.75.canonical_exons.bed")))
    setnames(exome_targets, c('chr', 'start', 'end', 'name', 'null', 'strand'))
    exome_targets$name = str_extract(exome_targets$name, '^[^:]+(?=:)')
    exome_targets[, null := NULL]
    exome_targets = exome_targets[chr %in% c(seq(1, 22), 'X')] # retain only autosomes and X
    exome_targets[, Hugo_Symbol := gsub(":.*$", "", name)]
    exome_targets[, Hugo_Symbol := ifelse(name %like% 'ENSG00000217792', 'LSP1P1', Hugo_Symbol)] # manual fix of bad annotation
    rm_genes = unique(exome_targets[, list(Hugo_Symbol, chr)])[, .N, keyby = Hugo_Symbol][N > 1] # these genes are mapping to multiple chromosomes
    exome_targets = exome_targets[Hugo_Symbol %ni% rm_genes$Hugo_Symbol]
    exome_targets
}

##

# Fetch a versioned freeze of OncoKB curated tumor-suppressors and oncogenes
oncokb = fromJSON(readLines(file.path(getSDIR(), 'data/OncoKB_allCuratedGenes_v1.19_patch_1.json'), warn=FALSE))
# If you have internet access, use the alternative below for the latest
## oncokb = fromJSON(readLines('http://oncokb.org/api/v1/genes', warn=FALSE))
oncokb_tsg = filter(oncokb, tsg=="TRUE") %>% select(hugoSymbol) %>% distinct(.)

#############################################
### definition of copy number calls in WGD
FACETS_CALL_table <- fread(file.path(getSDIR(), 'data/FACETS_CALL_table.tsv'))
setkey(FACETS_CALL_table, WGD, mcn, lcn)

## variable for filters, Shweta and Allison filters
fc_lu_table = FACETS_CALL_table %>% mutate(emtag = str_c(WGD,';',mcn,';',lcn))
## hey[, emtag := paste0(WGD,';',mcn,';',lcn)]
#############################################


annotate_integer_copy_number <- function(gene_level, amp_threshold){
    ### Portal/TCGA copy number calls as well as
    ### more advanced copy number calling, including CNLOH etc.

    ### lowest value of tcn for AMP
    AMP_thresh_tcn <- amp_threshold + 1

    ### assumes that WGD variable has already been set
    input_key <- key(gene_level)
    input_cols <- names(gene_level)
    gene_level[, mcn := tcn - lcn]
    setkey(gene_level, WGD, mcn, lcn)

    if("FACETS_CNA" %in% names(gene_level)) gene_level[, FACETS_CNA := NULL]
    if("FACETS_CALL" %in% names(gene_level)) gene_level[, FACETS_CALL := NULL]

    gene_level <- merge(gene_level, FACETS_CALL_table, sort = F, all.x = T)

    gene_level[tcn >= AMP_thresh_tcn, FACETS_CNA := 2]
    gene_level[tcn >= AMP_thresh_tcn, FACETS_CALL := "AMP"]

    output_cols <- c(input_cols, c("FACETS_CNA", "FACETS_CALL"))
    setcolorder(gene_level, output_cols)
    setkeyv(gene_level, input_key)
    gene_level
}


get_gene_level_calls <- function(cncf_files,
                                 gene_targets,
                                 WGD_threshold = 0.5, ### least value of frac_elev_major_cn for WGD
                                 amp_threshold = 5, ### total copy number greater than this value for an amplification
                                 max_seg_length = maxseg, ### genes in segments longer than this will be treated as diploid
                                 min_cf = min_cf_cutoff, ### genes in segments with cell fraction less than this will be treated as diploid
                                 mean_chrom_threshold = 0 ### total copy number also greater than this value multiplied by the chromosome mean for an amplification
                                 ){

    ### check version of data.table
    if(packageVersion("data.table") < "1.9.6"){stop("please update data.table to v1.9.6")}

    ### concatenate input files
    cncf_txt_list <- lapply(cncf_files, fread)
    names(cncf_txt_list) <- cncf_files
    ## concat_cncf_txt -- list of inputs as one data.table
    concat_cncf_txt <- rbindlist(cncf_txt_list, idcol = "filename", fill = TRUE)

    ### format concat_cncf_txt segment table
    concat_cncf_txt$chrom = as.character(concat_cncf_txt$chrom)
    concat_cncf_txt$Tumor_Sample_Barcode <- as.character(concat_cncf_txt$ID)
    concat_cncf_txt[, filename := NULL]
    concat_cncf_txt$chrom <- as.character(concat_cncf_txt$chrom)
    concat_cncf_txt[chrom == "23", chrom := "X"]
    setkey(concat_cncf_txt, chrom, loc.start, loc.end)

    if (!("tcn" %in% names(concat_cncf_txt))) {
        concat_cncf_txt[, c("tcn", "lcn") := list(tcn.em, lcn.em)]
    }

    ### estimate fraction of the genome with more than one copy from a parent
    ### a large fraction implies whole genome duplication
    concat_cncf_txt[, frac_elev_major_cn := sum(
        as.numeric((tcn - lcn) >= 2) * as.numeric(loc.end-loc.start), na.rm = T) /
      sum(as.numeric(loc.end-loc.start)
      ),
    by=Tumor_Sample_Barcode]

    ### Extract integer copy number for each probe from concat_cncf_txt
    fo_impact <- foverlaps(gene_targets, concat_cncf_txt, nomatch = NA,
                           by.x = c('chr', 'start', 'end'),
                           by.y = c('chrom', 'loc.start', 'loc.end'))
    fo_impact <- fo_impact[!is.na(ID)]
    fo_impact[,Hugo_Symbol:=name]

    if (!("cf" %in% names(fo_impact))) {
        fo_impact[, cf := cf.em]
    }

    ### Summarize copy number for each gene
    gene_level <- fo_impact[, list(chr = unique(chr),
                               seg.start=unique(loc.start),
                               seg.end=unique(loc.end),
                               frac_elev_major_cn=unique(frac_elev_major_cn),
                               Nprobes = .N),
                          keyby=list(Tumor_Sample_Barcode, Hugo_Symbol, tcn=tcn, lcn=lcn, cf=cf,
                                     tcn.em = tcn.em, lcn.em = lcn.em, cf.em = cf.em)]
#  }

### Collapsed into one call, see above
#   if(method == 'em'){

#       gene_level <- fo_impact[, list(chr = unique(chr),
#                              seg.start=unique(loc.start),
#                              seg.end=unique(loc.end),
# #                                start=unique(start),   ### with these uncommented, the per-gene summarization is broken (??)
# #                                end=unique(end),
#                              frac_elev_major_cn=unique(frac_elev_major_cn),
#                              Nprobes = .N),
#                         keyby=list(Tumor_Sample_Barcode, Hugo_Symbol, tcn=tcn.em, lcn=lcn.em)]
#   }

    ### fix bug where lcn == NA even when tcn is 1
    gene_level[tcn == 1 & is.na(lcn), lcn := 0]

    ### apply WGD threshold
    gene_level[, WGD := factor(ifelse(frac_elev_major_cn > WGD_threshold, "WGD", "no WGD"))]

    setkey(gene_level, Tumor_Sample_Barcode, Hugo_Symbol)

    ### focality requirement
    ### get (weighted) mean total copy number for the chromosome
#   mean_tcn_per_chr <- concat_cncf_txt[, list(mean_tcn =
#                                                sum(tcn * as.numeric(loc.end - loc.start)) /
#                                                sum(as.numeric(loc.end - loc.start))),
#                                       keyby=list(Tumor_Sample_Barcode, chrom)]
#   mean_tcn_per_chr_v <- with(mean_tcn_per_chr,
#                              structure(mean_tcn,
#                                        .Names = paste(Tumor_Sample_Barcode,
#                                                       chrom)))
#   gene_level[, mean_tcn_per_chr := mean_tcn_per_chr_v[paste(Tumor_Sample_Barcode, chr)]]

    ### annotate integer copy number
    gene_level <- annotate_integer_copy_number(gene_level, amp_threshold)

    ### remove duplicate entries for partial deletions
    ### the lower value is chosen on the basis that a
    ### partial deletion is a deletion but a
    ### partial amplification is not an amplification
    gene_level <- gene_level[order(FACETS_CNA, -Nprobes)]
    gene_level <- gene_level[!duplicated(gene_level, by=c("Tumor_Sample_Barcode", "Hugo_Symbol"))]

    ### Sort & set column classes
    gene_level[, FACETS_CNA := factor(FACETS_CNA, levels = c(-2:2))]
    gene_level[, chr := factor(chr, levels = c(1:22, "X", "Y"))]
    gene_level = as.data.table(gene_level)
    setkey(gene_level, Tumor_Sample_Barcode, Hugo_Symbol)


    ## Extract purity values
    ## March 2019, for filters
    ## within *_hisens.cncf.txt file, the purity is calculated as the maximum value not 1.0 and not NA
    ## For each unique Tumor_Sample_Barcode
    ##     subset cf.em; define purity as purity := max(cf.em)
    gene_level[cf.em < 1, purity:=max(cf.em), by=Tumor_Sample_Barcode]
    gene_level[, purity:=unique(purity)[!is.na(unique(purity))], by=Tumor_Sample_Barcode]

    ## Define CFcut, the cut-off for cell fraction, corrected for purity
    ## We suspect 60% of CF is pretty high, so that's the cutoff value set by min_cf_cutoff
    ## If purity is NA, then this cut-off will be ignored in downstream filters
    gene_level[, CFcut := min_cf * purity]

    ## Apply filters from Allison Richards and Shweta Chavan, March 2019
    ## Use FACETS EM as it gives more accurate results
    ##
    ## Decision tree:
    ## Inputs: AMP, AMP (LOH), HOMDEL
    ##
    ## For HOMDELs, keep if and only if size <10 Mb and # of genes <= 10
    ## Otherwise, throw awawy
    ## For AMPs and AMP (LOH) events, throw away if <10 Mb. Check TCN.em > 8, em.CF > 0.6 * purity. If both false, check # of gene <=10.
    ## If false, throw this event away.
    ##
    ## The 10 gene cut-off is based all annotated genes; for WES use
    ##

    ### Step 2.1 CONVERT GENELEVEL CALLS CNCF-based to EM-based
    genelevelcalls0 = gene_level  %>%
                     mutate(segid = str_c(Tumor_Sample_Barcode,';',chr,';',seg.start,';',seg.end)) %>%
                     mutate(mcn.em = tcn.em - lcn.em, seg.len = seg.end - seg.start)

    u.genelevelcalls = genelevelcalls0 %>%
                   select(Tumor_Sample_Barcode,tcn,lcn,tcn.em,lcn.em,chr,seg.start,seg.end,frac_elev_major_cn,WGD,mcn,FACETS_CNA,FACETS_CALL) %>% distinct(.) %>%
                   mutate(mcn.em = tcn.em - lcn.em, seg.len = seg.end - seg.start)

    frac_elev_major_cn.em = u.genelevelcalls %>%
                          group_by(Tumor_Sample_Barcode) %>%
                          summarize(frac_elev_major_cn.em = sum(as.numeric(mcn.em >= 2)*as.numeric(seg.len), na.rm = T) / sum(as.numeric(seg.len)))

    genelevelcalls0 = left_join(genelevelcalls0,frac_elev_major_cn.em, by="Tumor_Sample_Barcode")

    ### Step 2.2 Going back to original un-unique genelevel calls, just carry forward frac_elev_major_cn.em
    genelevelcalls0 = genelevelcalls0 %>%
                            mutate(WGD.em = ifelse(frac_elev_major_cn.em > WGD_threshold, "WGD", "no WGD")) %>%
                            mutate(emtag = str_c(WGD.em,';',mcn.em,';',lcn.em))

    genelevelcalls0 = genelevelcalls0 %>%
                            mutate(FACETS_CALL.em = plyr::mapvalues(emtag, fc_lu_table$emtag, fc_lu_table$FACETS_CALL)) %>% ##ASK
                            mutate(FACETS_CALL.em = ifelse(tcn.em >= 6,"AMP",FACETS_CALL.em)) %>%
                            mutate(FACETS_CALL.em = ifelse( (!is.na(tcn.em) & !is.na(lcn.em) & FACETS_CALL.em %ni% c(unique(fc_lu_table$FACETS_CALL)) ), "ILLOGICAL", FACETS_CALL.em))
    ##table(genelevelcalls0$FACETS_CALL.em)

    ##This line only works if processing exomes
    # 10 Genes cut-off is too large if processing IMPACT so need to count genes based on the exome bed file from ensembl, version 75
    # This section replicates what would happen if the whole exome bed was used instead of the IMPACT bed and counts the number of annotated genes in a given segment to be used for filtering in a later step
    if (args$targetFile == 'exome') {
        exome_bed = gene_targets
    } else {
        exome_bed = exome_targets()
    }
    cross <- foverlaps(exome_bed, concat_cncf_txt, nomatch = NA,
                       by.x = c('chr', 'start', 'end'),
                       by.y = c('chrom', 'loc.start', 'loc.end'))
    cross <- cross[!is.na(ID)]
    # cross[,Hugo_Symbol:=gsub(":.*$", "", name)]
    genecount <- cross[, list(chr = unique(chr),
                                 seg.start=unique(loc.start),
                                 seg.end=unique(loc.end),
                                 frac_elev_major_cn=unique(frac_elev_major_cn),
                                 Nprobes = .N),
                            keyby=list(Tumor_Sample_Barcode, Hugo_Symbol, tcn=tcn, lcn=lcn, cf=cf,
                                       tcn.em = tcn.em, lcn.em = lcn.em, cf.em = cf.em)]
    ### fix bug where lcn == NA even when tcn is 1
    genecount[tcn == 1 & is.na(lcn), lcn := 0]
    ### apply WGD threshold
    genecount[, WGD := factor(ifelse(frac_elev_major_cn > WGD_threshold, "WGD", "no WGD"))]
    setkey(genecount, Tumor_Sample_Barcode, Hugo_Symbol)
    ### annotate integer copy number
    genecount <- annotate_integer_copy_number(genecount, amp_threshold)
    ### remove duplicate entries for partial deletions
    ### the lower value is chosen on the basis that a
    ### partial deletion is a deletion but a
    ### partial amplification is not an amplification
    genecount <- genecount[order(FACETS_CNA, -Nprobes)]
    genecount <- genecount[!duplicated(genecount, by=c("Tumor_Sample_Barcode", "Hugo_Symbol"))]
    genecount = genecount  %>%
      mutate(segid = str_c(Tumor_Sample_Barcode,';',chr,';',seg.start,';',seg.end))
    seg.count = plyr::count(genecount, vars="segid") %>%
                            mutate(count = freq) %>% select(-c(freq))

    genelevelcalls0 = dplyr::inner_join(genelevelcalls0, seg.count, by = 'segid')

    ### Step 3 SET columns needed for filters & APPLY filters
    #genelevelcalls0= genelevelcalls0 %>% mutate(CFcut = plyr::mapvalues(Tumor_Sample_Barcode, gene_level$T, gene_level$CFcut))

    genelevelcalls0 = genelevelcalls0 %>%
        mutate(FACETS_CALL.ori = FACETS_CALL.em,
               FACETS_CALL.em = ifelse(FACETS_CALL.em %in% c("AMP","AMP (LOH)","AMP (BALANCED)","HOMDEL"),
                                        ifelse((FACETS_CALL.em %in% c("AMP","AMP (LOH)","AMP (BALANCED)") & seg.len < max_seg_length & (tcn.em > 8 | count <= 10 | ( !is.na(purity) & cf.em > CFcut ))), FACETS_CALL.em,
                                                ifelse( (FACETS_CALL.em == "HOMDEL" & seg.len < max_seg_length & count <= 10), FACETS_CALL.em, "ccs_filter")), FACETS_CALL.em ),
               review = FACETS_CALL.em == 'ccs_filter' & FACETS_CALL.ori == "HOMDEL" & Hugo_Symbol %in% unique(oncokb_tsg$hugoSymbol) & seg.len < 25000000,
               FACETS_CALL.em = FACETS_CALL.ori)
    ## table(genelevelcalls0$FACETS_CALL.em)

    # homdeltsg_review = filter(genelevelcalls0, FACETS_CALL.em == "ccs_filter", FACETS_CALL.ori == "HOMDEL", Hugo_Symbol %in% unique(oncokb_tsg$hugoSymbol), seg.len < 25000000)

    genelevelcalls0 = genelevelcalls0 %>%
        mutate(FACETS_CNA.em = plyr::mapvalues(FACETS_CALL.em, fc_lu_table$FACETS_CALL, fc_lu_table$FACETS_CNA)) %>%
        # mutate(FACETS_CNA.em = ifelse(FACETS_CALL.em=="ccs_filter",0,FACETS_CNA.em)) %>%
        mutate(FACETS_CNA.em = ifelse(FACETS_CALL.em=="ILLOGICAL",NA,FACETS_CNA.em)) %>%
        select(-c(FACETS_CNA, FACETS_CALL, CFcut, FACETS_CALL.ori, count, segid, seg.len, emtag, WGD.em, mcn, mcn.em, frac_elev_major_cn.em))
    
    genelevelcalls_final = genelevelcalls0 %>%
                            mutate(FACETS_CNA = FACETS_CNA.em, FACETS_CALL = FACETS_CALL.em) %>%
                            select(-c(FACETS_CALL.em,FACETS_CNA.em))


    # return(list(genelevelcalls_final = genelevelcalls_final, homdeltsg_review = homdeltsg_review)) ## return a list, and access these below
    filter(genelevelcalls_final, FACETS_CALL != 'DIPLOID')
}


convert_gene_level_calls_to_matrix_portal <- function(gene_level_calls){
    data_to_convert <- gene_level_calls[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "FACETS_CNA")]
    portal_output <- dcast(data_to_convert, Hugo_Symbol ~ Tumor_Sample_Barcode, value.var = "FACETS_CNA" )
    return(portal_output)
}

convert_gene_level_calls_to_matrix_ascna <- function(gene_level_calls){
    # gene_level_calls[is.na(gene_level_calls)] <- 'NA'
    gene_level_calls$ascna <- paste0(gene_level_calls$tcn,';',gene_level_calls$lcn)
    data_to_convert <- gene_level_calls[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "ascna")]
    ascna_output <- dcast(data_to_convert, Hugo_Symbol ~ Tumor_Sample_Barcode, value.var = "ascna" )
    ascna_output[is.na(ascna_output)] <- 'NA;NA'
    return(ascna_output)
}

#####################################################################################
#####################################################################################


if(!interactive()){

    parser = ArgumentParser()
    parser$add_argument('-f', '--filenames', type='character', nargs='+', help='list cncf files to be processed, to be concatenated into one R data.table.')
    parser$add_argument('-o', '--outfile', type='character', help='Output filename.')
    parser$add_argument('-t', '--targetFile', type='character', default='IMPACT468', help="IMPACT341/410/468, or a Picard interval list file of gene target coordinates [default IMPACT468]")
    parser$add_argument('-m', '--method', type='character', default='reg', help="If scna, creates a portal-friendly scna output file [default reg]")
    # parser$add_argument('-r', '--review_output_file', type='character', default='ccs_homdeltsg_review_candidates.txt', help="Output text file of canddiates for manual review")
    parser$add_argument('--min_cf_cutoff', type='double', default=0.6, help="The cell fraction cutoff such that genes in segments with cell fraction less than this will be treated as diploid")
    parser$add_argument('--max_seg_length', type='double', default=10000000, help="Genes in segments longer than this will be treated as diploid")
    args=parser$parse_args()

    filenames = args$filenames

    outfile = args$outfile
    method = args$method
    review_candidates = args$review_output_file
    min_cf_cutoff = args$min_cf_cutoff
    maxseg = args$max_seg_length

    ## Check that min_cf_cutoff is a fraction, i.e. 0<=min_cf_cutoff<=1, else throw an error
    if (min_cf_cutoff > 1 || min_cf_cutoff < 0){
        stop("min_cf_cutoff must be a fraction, a numeric between 0 and 1; please revise input min_cf_cutoff")
    }

    ## Check that max_seg_length is greater than 0
    if (maxseg < 0){
        stop("max_seg_length must be a numeric greater than 0; please revise input max_seg_length")
    }

    if (args$targetFile == "IMPACT341") {
        geneTargets = IMPACT341_targets()
    } else if (args$targetFile == "IMPACT410") {
        geneTargets = IMPACT410_targets()
    } else if (args$targetFile == "IMPACT468") {
        geneTargets = IMPACT468_targets()
    } else if (args$targetFile == 'exome') {
        geneTargets = exome_targets()
    } else {
        geneTargets <- suppressWarnings(fread(paste0('grep -v "^@" ',args$targetFile)))
        setnames(geneTargets, c("chr", "start", "end", "strand", "name"))
        setkey(geneTargets, chr, start, end)
        # If user gave us a DMP-IMPACT interval list, remove non-coding targets and reduce exon labels to gene names
        geneTargets <- geneTargets[grep("^Tiling_|^FP_|_intron_|_pseudo_|_MSI_|_tiling_|_promoter_", name, invert=TRUE),]
        geneTargets$name <- gsub("_target_.*$", "", geneTargets$name)
    }

    gene_level_calls = get_gene_level_calls(filenames, geneTargets)
    write.text(gene_level_calls, outfile)
    # This writes out homozygous deletions in known tumor suppressor genes that have been suppressed due to size or the number of genes
    # Known false negatives exist, such as RB1 that exists in a gene-rich region and so often (but not always) fails the 10 gene cut-off
    # This list is meant to be looked at manually to be sure that a true homdel is not suppressed (this will be improved in the next version)
    # fwrite(gene_level_calls$homdeltsg_review, paste0(gsub("[.]tsv|[.]txt","",outfile),"_TSG_ManualReview.txt"),row.names=FALSE, quote=FALSE, sep="\t")

    if(tolower(method) == 'scna'){
        scna_outfile = gsub(".txt", ".scna.txt", outfile)
        ascna_outfile = gsub(".txt", ".ascna.txt", outfile)
        portal_output = convert_gene_level_calls_to_matrix_portal(gene_level_calls)
        ascna_output = convert_gene_level_calls_to_matrix_ascna(gene_level_calls)
        write.text(portal_output, scna_outfile)
        write.text(ascna_output, ascna_outfile)
    }
}
