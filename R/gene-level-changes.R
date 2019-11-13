#' Get gene-level changes in copy number from FACETS output.
#'
#' This maps protein-coding genes onto copy-number segmentation from FACETS output derives the copy-number status from that, based on the integer copy number.
#' Additionally, the copy-number log-ratio of the SNPs falling insided the gene boundaries are used to perform a Z-test against the dipLogR baseline.
#' 
#' @param facets_output Full FACETS output from \code{run_facets}.
#' @param genome Genome build.
#' @param algorithm Choice between FACETS \code{em} and \code{cncf} algorithm.
#'
#' @return \code{data.frame} for all genes mapping onto a segment in the output segmentation, with the columns:
#' \itemize{
#'     \item{\code{median_cnlr_seg}:} {Median copy-number log-ratio for segment.}
#'     \item{\code{genes_on_seg}:} {Protein-coding genes on segment.}
#'     \item{\code{gene_snps}:} {Median copy-number log-ratio for segment.}
#'     \item{\code{gene_het_snps}:} {Median copy-number log-ratio for segment.}
#'     \item{\code{spans_segs}:} {Boolean indicating whether the gene boundaries map onto different segments.}
#'     \item{\code{cn_state}:} {Label for coyp-number configuration.}
#'     \item{\code{tsg}:} {Boolean indicating whether gene is a tumor suppressor.}
#' }
#' 
#' @import data.table
#' @importFrom plyr mapvalues
#' 
#' @examples
#' \dontrun{
#' gene_level_changes(test_facets_output, 'hg38', 'em')
#' }

#' @export
gene_level_changes = function(facets_output,
                              genome = c('hg19', 'hg38'),
                              algorithm = c('em', 'cncf')) {
    
    genome = match.arg(genome, c('hg19', 'hg38'), several.ok = FALSE)
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = FALSE)
    
    # Set variables
    segs = facets_output$segs
    snps = as.data.table(facets_output$snps)
    dipLogR = facets_output$dipLogR
    purity = facets_output$purity
    
    # Get WGD status
    fcna_output = calculate_fraction_cna(segs, ploidy, genome, algorithm)
    wgd = fcna_output$genome_doubled
    sample_ploidy = ifelse(wgd, round(ploidy), 2)
    
    # Map segments onto genes
    key_cols = c('chrom', 'start', 'end')
    genes = as.data.table(get(paste0('genes_', genome)))
    genes[, chrom := ifelse(chrom == 'X', 23, chrom)]
    genes = genes[chrom %in% seq(1, 23)][, chrom := as.integer(chrom)]
    
    setkeyv(genes, key_cols)
    
    segs = as.data.table(parse_segs(segs, algorithm))
    setkeyv(segs, key_cols)
    
    genes_segs = foverlaps(segs, genes, nomatch = 0, by.x = key_cols, by.y = key_cols) # i.start/i.end is the start/end interval of the segment
    
    # Count genes per segment
    genes_segs[, `:=` (genes_on_seg = .N), by = seg]
    
    # Map SNPs onto genes
    setDT(snps)[, `:=` (start = maploc, end = maploc)][, maploc := NULL]
    setkeyv(snps, key_cols)
    
    genes_snps = foverlaps(snps, genes, type = 'within', nomatch = 0, by.x = key_cols, by.y = key_cols)
    genes_snps = genes_snps[, list(
        # mean_cnlr = mean(cnlr),
        # sd_cnlr = sd(cnlr),
        # cnlr = list(cnlr),
        snps = .N,
        het_snps = sum(het == 1),
        seg = seg[which.max(start-end)], # this selects the segment which represents the larger region
        spans_segs = length(unique(seg)) > 1
    ), keyby = gene]
    
    # Combine gene-segmentation and gene-SNP mapping
    # Select segment from the SNP mapping
    genes_all = merge(genes_segs, genes_snps, by.x = 'gene', by.y = 'gene')[seg.x == seg.y]
        
    # Map to copy-number states
    genes_all[, cn_state := mapvalues(paste(wgd, tcn-lcn, lcn, sep = ':'),
                                      copy_number_states$map_string, copy_number_states$call,
                                      warn_missing = FALSE)]
    genes_all[, cn_state := ifelse(!cn_state %in% copy_number_states$call, 'INDETERMINATE', cn_state)]
    
    # Test on cnlr against baseline
    # cn0_dipLogR = unique(segs$cnlr.median.clust)[which.min(abs(unique(segs$cnlr.median.clust)-dipLogR))]
    # cn0_segs = segs[cnlr.median.clust == cn0_dipLogR]
    # cn0_snps = snps[segclust %in% cn0_segs$segclust]
    # cn0_snps = cn0_snps[between(cnlr, quantile(cn0_snps$cnlr, .25), quantile(cn0_snps$cnlr, .75))] # remove noise
    
    # Perform Z test, calculate fold change
    # genes_all[, mean_cnlr := mean(unlist(cnlr)), by = seq_len(nrow(genes_all))][, `:=` (
    #     pval = ifelse(mean_cnlr > mean(cn0_snps$cnlr),
    #                   two_sample_z(cn0_snps$cnlr, unlist(cnlr)),
    #                   two_sample_z(unlist(cnlr), cn0_snps$cnlr)),
    #     fold_change = ifelse(mean_cnlr - mean(cn0_snps$cnlr) < 0,
    #                          -2^(-(mean_cnlr - mean(cn0_snps$cnlr))),
    #                          2^(mean_cnlr - mean(cn0_snps$cnlr)))
    # ), by = seq_len(nrow(genes_all))]
    
    # Add filter flags
    max_gene_count = 10
    max_seg_length = 1e7
    min_ccf = 0.6 * purity
    genes_all[, filter := 'PASS']
    genes_all[cn_state %like% '^AMP' & length > max_seg_length & (tcn > 8 | genes_on_seg > max_gene_count) , filter := add_tag(filter, 'unfocal_amp')]
    genes_all[cn_state %like% '^AMP' & length > max_seg_length & (!is.na(purity) & cf < min_ccf) , filter := add_tag(filter, 'subclonal_amp')]
    genes_all[cn_state == 'HOMDEL' & length > max_seg_length & genes_on_seg > max_gene_count, filter := add_tag(filter, 'unfocal_del')]
    genes_all[cn_state == 'HOMDEL' & tsg == TRUE & length < max_seg_length, filter := 'RESCUE']
    
    # Clean up data frame
    genes_all = genes_all[, setnames(.SD,
                                     c('seg.x', 'start', 'end', 'i.start', 'i.end', 'length', 'snps', 'het_snps', 'cnlr.median'),
                                     c('seg', 'gene_start', 'gene_end', 'seg_start', 'seg_end', 'seg_length', 'gene_snps', 'gene_het_snps', 'median_cnlr_seg'))]
    genes_all[, `:=` (num.mark = NULL, nhet = NULL, mafR = NULL, mafR.clust = NULL, seg.y = NULL, cnlr.median.clust = NULL)]
    
    data.frame(genes_all)
}

# Help functions --------------------------------------------------------------------------------------------------

# two_sample_z = function(a, b) {
#     if (length(a) < 5 | length(b) < 5) {
#         NA_real_
#     } else {
#         se_a = sd(a)/sqrt(length(a))
#         se_b = sd(b)/sqrt(length(b))
#         se = sqrt(se_a^2 + se_b^2)
#         z = (mean(a)-mean(b))/se
#         pnorm(z, lower.tail = TRUE)
#     }
# }

add_tag = function(filter, tag) {
    ifelse(filter == 'PASS',
           tag,
           paste(filter, tag, sep = ';'))
}


