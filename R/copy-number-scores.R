#' Calculate copy-number based scores
#'
#' @param segs
#' @param ploidy
#' @param genome
#' @param algorithm
#'
#' @return
#' @name copy_number_scores
NULL
#' 
#' @export
#' @rdname copy_number_scores
calculate_fcna = function(segs,
                          ploidy,
                          genome = c('hg18', 'hg19', 'grch38'),
                          algorithm = c('em', 'cncf')) {

    genome = match.arg(genome)
    algorithm = match.arg(algorithm)

    # Centromere locations
    genome = get(genome)

    segs = parse_segs(segs, algorithm)

    # Create chrom_info for sample
    sample_chrom_info = get_sample_genome(segs, genome)

    # Calculated length of interrogated genome
    interrogated_genome = sum(sample_chrom_info$size)
    autosomal_genome = sum(sample_chrom_info$size[sample_chrom_info$chr %in% 1:22])

    # Check for whole-genome duplication // PMID 30013179
    wgd_treshold = 0.5 # treshold
    frac_elevated_mcn = sum(segs$length[which(segs$mcn >= 2 & segs$chrom %in% 1:22)])/autosomal_genome
    wgd = frac_elevated_mcn > wgd_treshold

    # Calculate fraction of genome altered
    sample_ploidy = ifelse(wgd, round(ploidy), 2)
    if (!wgd) {
        diploid_length = sum(segs$length[which(segs$tcn == sample_ploidy & segs$lcn == 1)])
    } else if (wgd) {
        diploid_length = sum(segs$length[which(segs$tcn == sample_ploidy & segs$lcn >= 1)])
    }
    frac_altered = (interrogated_genome-diploid_length)/interrogated_genome

    # Altered arms, slightly adjusted to PMID 29622463
    segs = left_join(segs, select(sample_chrom_info, chr, centromere), by = c('chrom' = 'chr'))
    altered_arms = sapply(unique(segs$chrom), function(x) {
        seg_p = segs[which(segs$chrom == x & segs$start < segs$centromere),]
        seg_p$end = ifelse(seg_p$end > seg_p$centromere, seg_p$centromere, seg_p$end)
        seg_p$length = seg_p$end - seg_p$start
        seg_q = segs[which(segs$chrom == x & segs$end > segs$centromere),]
        seg_q$start = ifelse(seg_q$start < seg_q$centromere, seg_q$centromere, seg_q$start)
        seg_q$length = seg_q$end - seg_q$start

        if (!wgd) {
            seg_p_unaltered = sum(seg_p$length[which(seg_p$tcn == sample_ploidy & seg_p$lcn == 1)])
            seg_q_unaltered = sum(seg_q$length[which(seg_q$tcn == sample_ploidy & seg_q$lcn == 1)])
        } else if (wgd) {
            seg_p_unaltered = sum(seg_p$length[which(seg_p$tcn == sample_ploidy & seg_p$lcn >= 1)])
            seg_q_unaltered = sum(seg_q$length[which(seg_q$tcn == sample_ploidy & seg_q$lcn >= 1)])
        }
        paste0(paste0(x, c('p', 'q'))[c((sample_chrom_info$plength[sample_chrom_info$chr == x]-seg_p_unaltered)/sample_chrom_info$plength[sample_chrom_info$chr == x] > .8,
                                        (sample_chrom_info$qlength[sample_chrom_info$chr == x] - seg_q_unaltered)/sample_chrom_info$qlength[sample_chrom_info$chr == x] > .8) %>%
                                          map_lgl(~ifelse(is.na(.), FALSE, .))], collapse = ',')
    }) %>% strsplit(., ',') %>%
        unlist %>%
        discard(. %in% c('', '13p', '14p', '15p', '21p', '22p'))

    # weighted fraction copy-number altered
    frac_altered_w = select(sample_chrom_info, chr, p = plength, q = qlength) %>%
        gather(arm, length, -chr) %>%
        filter(!paste0(chr, arm) %in% c('13p', '14p', '15p', '21p', '22p')) %>%
        summarize(sum(length[paste0(chr, arm) %in% altered_arms])/sum(length)) %>%
        as.numeric()

    data.frame(
        wgd = wgd,
        fcna = frac_altered,
        weighted_fcna = frac_altered_w,
        aneuploidy_score = length(altered_arms),
        altered_arms = paste0(altered_arms, collapse = ',')
    )
}

parse_segs = function(segs, algorithm) {
    if (algorithm == 'em') {
        segs$tcn = segs$tcn.em
        segs$lcn = segs$lcn.em
    }

    mutate(segs,
           length = end - start,
           lcn = ifelse(tcn <= 1, 0, lcn),  # correct mcn, lcn for cases of tcn = 1 // sometimes FACETS set lcn = NA when tcn = 1, when it clearly has to be 0
           mcn = tcn - lcn)
}

get_sample_genome = function(segs, genome) {
    ldply(unique(segs$chrom), function(x) {
        centromere = genome$centromere[which(genome$chrom == x)]
        chrstart = min(segs$start[which(segs$chrom == x)])
        chrend = max(segs$end[which(segs$chrom == x)])
        size = chrend-chrstart
        plength = centromere - chrstart
        if (plength < 0) plength = 0 # acrocentric chromosomes
        qlength = chrend - centromere
        c(chr = x,
          centromere = centromere,
          chrstart = chrstart,
          chrend = chrend,
          size = size,
          plength = plength,
          qlength = qlength)
    })
}
