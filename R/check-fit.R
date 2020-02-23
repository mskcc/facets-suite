#' Sample QC
#' 
#' Generate QC metrics for sample. Can use mutation calls in MAF file.
#'
#' @param facets_output Output from \code{run_facets}.
#' @param genome Genome build.
#' @param algorithm Choose assessing the fit from the \code{em} or \code{cncf} algorithm.
#' @param maf Optional: mutation calls for assessed samples, should only include mutations for given sample.
#'
#' @return A list object with the following items:
#' \itemize{
#'       \item{\code{dipLogR_flag}:} {Boolean indicating extreme dipLogR value.}
#'       \item{\code{n_alternative_dipLogR}:} {Number of alternative dipLogR values.}
#'       \item{\code{n_dip_bal_segs}, \code{frac_dip_bal_segs}:} {Number of balanced segments at dipLogR and the fraction of genome they represent.}
#'       \item{\code{n_dip_imbal_segs}, \code{frac_dip_imbal_segs}:} {Number of imbalanced segments at dipLogR and the fraction of genome they represent.}
#'       \item{\code{n_amp}:} {Number of segments at total copy number >= 10.}
#'       \item{\code{n_homdels}:} {Number of homozygously deleted segments (total copy number = 0).}
#'       \item{\code{n_homdels_clonal}, \code{frac_homdels_clonal}:} {Number of clonal homdel segments and the fraction of the genome they represent.}
#'       \item{\code{n_cn_states}:} {Number of unique copy-number states (i.e. combinations of major and minor copy number).}
#'       \item{\code{n_segs}:} {Number of segments.}
#'       \item{\code{n_cnlr_clusters}:} {Number of copy-number log-ratio clusters}
#'       \item{\code{n_lcn_na}:} {Number of segments where no minor copy number was inferred (lcn is NA).}
#'       \item{\code{n_loh}, \code{n_loh}:} {Number of segments where the minor copy number is 0 and the fraction of the genome they represent.}
#'       \item{\code{n_snps}:} {Number of SNPs used for segmentation.}
#'       \item{\code{n_het_snps}, \code{frac_het_snps}:} {Number of heterozyous SNPs used for segmentation and their fraction of the total.}
#'       \item{\code{n_het_snps_hom_in_tumor_1pct}, \code{frac_het_snps_hom_in_tumor_1pct}:} {Number of heterozyous SNPs where the tumor allele frequency is <0.01/>0.99 their fraction of the total.}
#'       \item{\code{n_het_snps_hom_in_tumor_5pct}, \code{frac_het_snps_hom_in_tumor_5pct}:} {Number of heterozyous SNPs where the tumor allele frequency is <0.05/>0.95 their fraction of the total.}
#'       \item{\code{mean_cnlr_residual}, \code{sd_cnlr_residual}:} {Mean and standard deviation of SNPs' log-ratio from their segments copy-number log-ratio.}
#'       \item{\code{n_segs_discordant_tcn}, \code{frac_segs_discordant_tcn}:} {Number of segments where the naïve and EM algorithm estimates of the total copy number are discordant and the fraction of the genome they represent.}
#'       \item{\code{n_segs_discordant_lcn}, \code{frac_segs_discordant_lcn}:} {Number of segments where the naïve and EM algorithm estimates of the minor copy number are discordant and the fraction of the genome they represent.}
#'       \item{\code{n_segs_discordant_both}, \code{frac_segs_discordant_both}:} {Number of segments where the naïve and EM algorithm estimates of the both copy numbers are discordant and the fraction of the genome they represent.}
#'       \item{\code{n_segs_icn_cnlor_discordant}, \code{frac_icn_cnlor_discordant}:} {Number of clonal segments where the log-ratio shows balance but the copy-number solution does not, and the reverse, and the fraction of the genome they represent.}
#'       \item{\code{dip_median_vaf}:} {If MAF input: median tumor VAF of somatic mutations on clonal segments with total copy number 2 and allelic balance.}
#'       \item{\code{n_homdel_muts}:} {If MAF input: number of somatic mutations in homozygously deleted segments.}
#'       \item{\code{median_vaf_homdel_muts}:} {If MAF input: Median tumor VAF of somatic  mutations homozygously deleted segments.}
#'    }
#'
#' @importFrom dplyr distinct
#' @importFrom purrr map_if
#' @import data.table

#' @export
check_fit = function(facets_output,
                     genome = c('hg19', 'hg18', 'hg38'),
                     algorithm = c('em', 'cncf'),
                     maf = NULL) {
    
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = F)
    genome_choice = get(match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = F))
    
    # Set variables
    segs = as.data.table(facets_output$segs)
    snps = facets_output$snps
    dipLogR = facets_output$dipLogR
    alballogr = as.numeric(facets_output$alballogr[, 1])
    purity = facets_output$purity
    
    fcna_output = calculate_fraction_cna(segs, facets_output$ploidy, genome, algorithm)
    wgd = fcna_output$genome_doubled
    
    segs[, `:=` (tcn.original = tcn, lcn.original = lcn, cf.original = cf,
                 tcn.em.original = tcn.em, lcn.em.original = lcn.em, cf.em.original = cf.em)] # these are retained since the parse_segs function consolidates cncf/em solutions
    segs = parse_segs(segs, algorithm)
    
    segs = as.data.table(segs)
    setkey(segs, chrom, start, end)
    snps = as.data.table(snps)
    setkey(snps, chrom, maploc)
    
    # Label clonal segments
    segs[, clonal := cf >= (purity * 0.8)]
    
    # Subset on autosomes
    auto_segs = segs[chrom < 23]
    
    # Check for extrema dipLogR values
    dipLogR_flag = abs(facets_output$dipLogR) > 1
    
    # Check for alternative dipLogR values
    n_alt_dipLogR = length(setdiff(alballogr, dipLogR))
    
    # Check for balance/imbalance of copy-number at segments at dipLogR
    # cnlr.median.clust == dipLogR for purity runs, for others, find the closest one
    cnlr_clusts = unique(segs$cnlr.median.clust)
    cnlr_clust_value = cnlr_clusts[which.min(abs(cnlr_clusts - dipLogR))]
    
    dip_bal_segs = auto_segs[cnlr.median.clust == cnlr_clust_value & mcn == lcn, ]
    n_dip_bal_segs = nrow(dip_bal_segs)
    frac_dip_bal_segs = sum(dip_bal_segs$length)/sum(auto_segs$length)
    
    dip_imbal_segs = auto_segs[cnlr.median.clust == cnlr_clust_value & mcn != lcn & !is.na(lcn), ]
    n_dip_imbal_segs = nrow(dip_imbal_segs)
    frac_dip_imbal_segs = sum(dip_imbal_segs$length)/sum(auto_segs$length)
    
    #############################
    ### Other metrics
    #############################
    ## fraction genome unaltered (2-1) -- too high --> low purity?
    segs_unaltered = auto_segs[tcn == 2 & lcn == 1, ]    
    n_segs_unaltered = nrow(segs_unaltered)    
    frac_genome_unaltered = sum(segs_unaltered$length)/sum(auto_segs$length)
    
    segs_below_dipLogR = auto_segs[cnlr.median.clust < dipLogR]    # part of ploidy filter
    n_segs_below_dipLogR = nrow(segs_below_dipLogR)    # part of ploidy filter
    frac_below_dipLogR = sum(segs_below_dipLogR$length)/sum(auto_segs$length)
    
    segs_balanced_odd_tcn = auto_segs[mafR <= 0.025 & (tcn %% 2) != 0 & clonal ==T & tcn < 8,]
    n_segs_balanced_odd_tcn = nrow(segs_balanced_odd_tcn)
    frac_balanced_odd_tcn = sum(segs_balanced_odd_tcn$length)/sum(auto_segs$length)
    
    segs_imbalanced_diploid_cn = auto_segs[clonal==T & mafR > 0.1 & !is.na(lcn) & lcn != 0  & (as.double(tcn) / lcn) == 2, ]
    n_segs_imbalanced_diploid_cn = nrow(segs_imbalanced_diploid_cn)
    frac_imbalanced_diploid_cn = sum(segs_imbalanced_diploid_cn$length)/sum(auto_segs$length)
    
    segs_lcn_greater_mcn = auto_segs[clonal==T & lcn < mcn,]
    n_segs_lcn_greater_mcn = nrow(segs_lcn_greater_mcn)
    frac_lcn_greater_mcn = sum(segs_lcn_greater_mcn$length)/sum(auto_segs$length)
    
    mafr_median_all = median(auto_segs$mafR, na.rm = T)
    mafr_median_clonal = median(auto_segs[clonal==T, ]$mafR, na.rm = T)
    mafr_n_gt_1 = nrow(auto_segs[mafR > 1, ])
    
    # Number of high-level amplifications and homozygous deletions
    # Clonal homdels, how much of the genome do they represent
    n_amps = nrow(segs[tcn >= 10])
    homdels = auto_segs[tcn == 0]
    n_homdels = nrow(homdels)
    
    clonal_homdels = auto_segs[tcn == 0 & clonal == TRUE]
    n_homdels_clonal = nrow(clonal_homdels)
    
    frac_homdels = sum(homdels$length)/sum(auto_segs$length)
    frac_homdels_clonal = sum(clonal_homdels$length)/sum(auto_segs$length)
    
    # Count number of unique copy-number states and total number of segments
    n_cn_states = nrow(distinct(segs, tcn, lcn))
    n_segs = nrow(segs)
    n_cnlr_clusters = length(unique(segs$cnlr.median.clust))
    
    # Number of segments with lcn==NA
    n_lcn_na = nrow(auto_segs[is.na(lcn)])
    frac_lcn_na = sum(auto_segs[is.na(lcn)]$length)/sum(auto_segs$length)
    
    # Number of segments with LOH
    loh_segs = auto_segs[lcn == 0,]
    n_loh = nrow(loh_segs)
    frac_loh = sum(loh_segs$length)/sum(auto_segs$length)
    
    # Check fraction of subclonal events within autosomal chromosomes
    subclonal_segs = auto_segs[clonal == F, ]
    n_segs_subclonal = nrow(subclonal_segs)
    frac_segs_subclonal = sum(subclonal_segs$length)/sum(auto_segs$length)
    
    # Compile SNP statistics
    n_snps = nrow(snps)
    het_snps = snps[het == 1]
    n_het_snps = nrow(het_snps)
    frac_het_snps = n_het_snps/n_snps
    
    n_snps_with_300x_in_tumor = nrow(snps[rCountT > 300, ])
    n_het_snps_with_300x_in_tumor = nrow(het_snps[rCountT > 300, ])
    n_het_snps_hom_in_tumor_1pct = nrow(het_snps[  (vafT < 0.01 | vafT > 0.99), ])
    n_het_snps_hom_in_tumor_5pct = nrow(het_snps[ (rCountN > 35 & rCountT > 35) & (vafT < 0.05 | vafT > 0.95), ])
    frac_het_snps_hom_in_tumor_1pct = n_het_snps_hom_in_tumor_1pct/n_het_snps
    frac_het_snps_hom_in_tumor_5pct = n_het_snps_hom_in_tumor_5pct/n_het_snps
    
    # Check the mean/standard deviation of the cnlr
    snps[, cnlr_residual := cnlr - median(cnlr), by = seg]
    mean_cnlr_residual = mean(snps$cnlr_residual)
    sd_cnlr_residual = sd(snps$cnlr_residual)
    
    # Check concordance between CNCF and EM fits 
    discordant_segs = auto_segs[, `:=` (
        discordant_tcn = (tcn.em.original != tcn.original | (is.na(tcn.em.original) | is.na(tcn.original))) & !(is.na(tcn.em.original) & is.na(tcn.original)),
        discordant_lcn = (lcn.em.original != lcn.original | (is.na(lcn.em.original) | is.na(lcn.original))) & !(is.na(lcn.em.original) & is.na(lcn.original))
    )][(discordant_tcn == TRUE | discordant_lcn == TRUE) & tcn < 10]
    discordant_stats = discordant_segs[, list(
        n_discordant_tcn = sum(discordant_tcn & !discordant_lcn),
        length_discordant_tcn = sum(length[discordant_tcn & !discordant_lcn]),
        n_discordant_lcn = sum(length[!discordant_tcn & discordant_lcn]),
        length_discordant_lcn = sum(length[!discordant_tcn & discordant_lcn]),
        n_discordant_both = sum(discordant_tcn & discordant_lcn),
        length_discordant_both = sum(length[discordant_tcn & discordant_lcn])
    )]
    evaluable_length = segs[chrom <= 22 & tcn < 10][, list(
        tcn = sum(length[!(is.na(tcn.em.original) & is.na(tcn.original))]),
        lcn = sum(length[!(is.na(lcn.em.original) & is.na(lcn.original))]),
        both = sum(length[!(is.na(tcn.em.original) & is.na(tcn.original)) & !(is.na(lcn.em.original) & is.na(lcn.original))])
    )]
    
    # Check concordance between log odds-ratio and integer copy-number with regards to balance
    # Count segments where logOR shows balanced but tcn/lcn is imbalanced and vice-versa.
    clonal_segs = auto_segs[clonal == TRUE & !is.na(lcn)]
    clonal_segs_disc_icn = clonal_segs[, `:=` (
        icn_bal_mafr_high = lcn == mcn & mafR > 0.05,
        icn_imbal_mafr_low = lcn != mcn & mafR < 0.05
    )][icn_bal_mafr_high == TRUE | icn_imbal_mafr_low == TRUE]
    
    n_icn_cnlor_discordant = nrow(clonal_segs_disc_icn)
    frac_icn_cnlor_discordant = sum(clonal_segs_disc_icn$length)/sum(auto_segs$length)
    
    # Output all values
    # n: denotes number
    # frac: denotes fraction of assesses genome
    output = list(
        dipLogR_flag = dipLogR_flag,
        n_alternative_dipLogR = n_alt_dipLogR,
        wgd = wgd,
        n_dip_bal_segs = n_dip_bal_segs,
        frac_dip_bal_segs = frac_dip_bal_segs,
        n_dip_imbal_segs = n_dip_imbal_segs,
        frac_dip_imbal_segs = frac_dip_imbal_segs,
        n_amps = n_amps,
        n_homdels = n_homdels,
        frac_homdels = frac_homdels,
        n_homdels_clonal = n_homdels_clonal,
        frac_homdels_clonal = frac_homdels_clonal,
        n_cn_states = n_cn_states,
        n_segs = n_segs,
        n_cnlr_clusters = n_cnlr_clusters,
        n_lcn_na = n_lcn_na,
        n_loh = n_loh,
        frac_loh = frac_loh,
        n_segs_subclonal = n_segs_subclonal,
        frac_segs_subclonal = frac_segs_subclonal,
        n_segs_below_dipLogR = n_segs_below_dipLogR,
        frac_below_dipLogR = frac_below_dipLogR,
        n_segs_balanced_odd_tcn = n_segs_balanced_odd_tcn,
        frac_balanced_odd_tcn = frac_balanced_odd_tcn,
        n_segs_imbalanced_diploid_cn = n_segs_imbalanced_diploid_cn,
        frac_imbalanced_diploid_cn = frac_imbalanced_diploid_cn,
        n_segs_lcn_greater_mcn = n_segs_lcn_greater_mcn,
        frac_lcn_greater_mcn = frac_lcn_greater_mcn,
        n_snps = n_snps,
        n_het_snps = n_het_snps,
        frac_het_snps = frac_het_snps,
        n_snps_with_300x_in_tumor = n_snps_with_300x_in_tumor,
        n_het_snps_with_300x_in_tumor = n_het_snps_with_300x_in_tumor,
        n_het_snps_hom_in_tumor_1pct = n_het_snps_hom_in_tumor_1pct,
        n_het_snps_hom_in_tumor_5pct = n_het_snps_hom_in_tumor_5pct,
        frac_het_snps_hom_in_tumor_1pct = frac_het_snps_hom_in_tumor_1pct,
        frac_het_snps_hom_in_tumor_5pct = frac_het_snps_hom_in_tumor_5pct,
        mean_cnlr_residual = mean_cnlr_residual,
        sd_cnlr_residual = sd_cnlr_residual,
        n_segs_discordant_tcn = discordant_stats$n_discordant_tcn,
        frac_discordant_tcn = discordant_stats$length_discordant_tcn/evaluable_length$tcn,
        n_segs_discordant_lcn = discordant_stats$n_discordant_lcn,
        frac_discordant_lcn = discordant_stats$length_discordant_lcn/evaluable_length$lcn,
        n_segs_discordant_both = discordant_stats$n_discordant_both,
        frac_discordant_both = discordant_stats$length_discordant_both/evaluable_length$both,
        n_segs_icn_cnlor_discordant = n_icn_cnlor_discordant,
        frac_icn_cnlor_discordant = frac_icn_cnlor_discordant,
        mafr_median_all = mafr_median_all,
        mafr_median_clonal = mafr_median_clonal,
        mafr_n_gt_1 = mafr_n_gt_1
    )
    
    
    
    # If input MAF is provided add some stats based on this
    if (!is.null(maf)) {
        maf = as.data.table(maf)
        
        maf[, `:=` (
            Chromosome = ifelse(Chromosome == 'X', 23, Chromosome),
            t_var_freq = t_alt_count/(t_alt_count+t_ref_count)
        )][, Chromosome := as.integer(Chromosome)]
        setkey(maf, Chromosome, Start_Position, End_Position)
        maf = foverlaps(maf, segs, mult = 'first', nomatch = NA,
                        by.x = c('Chromosome', 'Start_Position', 'End_Position'),
                        by.y = c('chrom', 'start', 'end'))
        
        # Median mutation VAF at clonal 2:1 segments
        output$dip_median_vaf = maf[clonal == TRUE & mcn == 1 & lcn == 1][, median(t_var_freq)]
        
        # Mutations at homdels
        homdel_muts = maf[tcn == 0]
        output$n_homdel_muts = nrow(homdel_muts)
        output$median_vaf_homdel_muts = median(homdel_muts$t_var_freq)
    }
    
    # Return rounded values for fractions
    map_if(output, is.double, .f = function(x) signif(x, 2))
}

