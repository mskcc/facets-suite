#' Copy-number based scores
#' 
#' Calculate the following 
#' \describe{
#'  \item{\strong{Fraction of genome altered:}}{Fraction of genome altered and genome doubling flag.}
#'  \item{\strong{Fraction LOH:}}{Fraction of genome with LOH and flag for hypoploidy.}
#'  \item{\strong{LST score:}}{Large-scale state transitions, see source URL.}
#'  \item{\strong{NtAI:}}{Telomeric allelic imbalance, see source URL.}
#'  \item{\strong{HRD-LOH:}}{HRD-LOH score, see source URL.}
#' }
#'
#' @param segs FACETS segmentation output.
#' @param snps FACETS SNP output.
#' @param ploidy Sample ploidy.
#' @param genome Genome build.
#' @param algorithm Choice between FACETS \code{em} and \code{cncf} algorithm.
#' @param min_size Minimum length of segment, as defined per function.
#' @param min_probes Minimum number of SNPs per segment, as defined per function.
#' @param hypoploidy_threshold Threshold for hypoploid call.
#' 
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/26015868}
#'
#' @return List with one or more values from function.
#'
#' @importFrom dplyr left_join summarize filter mutate group_by pull %>% select matches ends_with
#' @importFrom purrr map_dfr
#' @importFrom diptest dip.test
#' 
#' @examples 
#' calculate_fraction_cna(test_facets_output$segs, test_facets_output$ploidy, 'hg38', 'em')
#' calculate_loh(test_facets_output$segs, test_facets_output$snps, 'hg38', 'em')
#' calculate_lst(test_facets_output$segs, test_facets_output$ploidy, 'hg38', 'em')
#' calculate_ntai(test_facets_output$segs, test_facets_output$ploidy, 'hg38', 'em')
#' calculate_hrdloh(test_facets_output$segs, test_facets_output$ploidy, 'em')
#' 
#' @name copy_number_scores
NULL

#' @export 
#' @rdname copy_number_scores
calculate_fraction_cna = function(segs,
                                  ploidy,
                                  genome = c('hg19', 'hg18', 'hg38'),
                                  algorithm = c('em', 'cncf')) {
    
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = FALSE)
    
    # Centromere locations
    genome = get(match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE))
    
    # Create chrom_info for sample
    segs = parse_segs(segs, algorithm)
    sample_chrom_info = get_sample_genome(segs, genome)
    
    wgd = is_genome_doubled(segs, sample_chrom_info, treshold = 0.5)
    
    # Calculate fraction of genome altered
    interrogated_genome = sum(as.numeric(sample_chrom_info$size))
    if (!wgd) {
        diploid_length = sum(as.numeric(segs$length[which(segs$tcn == 2 & segs$lcn == 1)]))
    } else if (wgd) {
        diploid_length = sum(as.numeric(segs$length[which((segs$tcn == 2 & segs$lcn == 1) | (segs$tcn == 4 & segs$lcn == 2))]))
    }
    frac_altered = (interrogated_genome - diploid_length) / interrogated_genome
    
    list(
        genome_doubled = wgd,
        fraction_cna = frac_altered
    )
}

#' @export 
#' @rdname copy_number_scores
calculate_loh = function(segs,
                         snps,
                         genome = c('hg19', 'hg18', 'hg38'),
                         algorithm = c('em', 'cncf'),
                         hypoploidy_threshold = 0.5) {
    
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = FALSE)
    
    # Centromere locations
    genome = get(match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE))
    
    # Create chrom_info for sample
    segs = parse_segs(segs, algorithm)
    segs = filter(segs, chrom %in% 1:22)
    
    # Fraction LOH
    frac_loh = group_by(segs, chrom) %>% 
        summarize(chrom_length = as.numeric(sum(length[!is.na(lcn)])),
                  loh = as.numeric(sum(length[lcn == 0]))) %>% 
        mutate(loh = ifelse(is.na(loh), 0, loh)) %>% 
        summarize(fraction_loh = sum(loh, na.rm = TRUE)/sum(chrom_length),
                  loh_chromosomes = paste(c(chrom[(loh / chrom_length) > 0.75]), collapse = ','))
    
    # Check for true hypoploidy
    het_snps = filter(snps, het == 1) %>% 
        left_join(segs, ., by = 'seg')
    loh_snps = het_snps[het_snps$lcn == 0, ]
    bal_snps = het_snps[which(het_snps$lcn == (het_snps$tcn / 2)), ]
    gain_snps = het_snps[which(het_snps$lcn != (het_snps$tcn / 2) & het_snps$lcn > 0), ]
    loh_test = diptest::dip.test(loh_snps$valor)$p.value
    bal_test = diptest::dip.test(bal_snps$valor)$p.value
    gain_test = diptest::dip.test(gain_snps$valor)$p.value
    
    if (loh_test < .05 & bal_test > 0.05 &
        length(loh_snps) > 10 & length(bal_snps) > 10) {
        hypo_test = TRUE
    } else if (loh_test < .05 & bal_test > 0.05 &
               length(loh_snps) > 10 & length(bal_snps) > 10) {
        hypo_test = TRUE
    } else {
        hypo_test = F
    }
    
    list(
        fraction_loh = frac_loh$fraction_loh,
        loh_chromosomes = frac_loh$loh_chromosomes,
        hypoploid = hypo_test == TRUE & frac_loh$fraction_loh >= hypoploidy_threshold
    )
}

#' @export 
#' @rdname copy_number_scores
calculate_ntai = function(segs,
                          ploidy,
                          genome = c('hg19', 'hg18', 'hg38'),
                          algorithm = c('em', 'cncf'),
                          min_size = 0,
                          min_probes = 250) {
    
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = FALSE)
    
    # Centromere locations
    genome = get(match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE))
    
    # Create chrom_info for sample
    segs = parse_segs(segs, algorithm) %>% 
        filter(chrom %in% 1:22)
    sample_chrom_info = get_sample_genome(segs, genome)
    
    # Filter segments that do not pass filters
    segs = segs[which(segs$num.mark >= min_probes & segs$length > min_size), ]
    
    # Shrink segments with identical copy number
    segs = map_dfr(unique(segs$chrom), function(i) {
        chrom_seg = segs[which(segs$chrom == i), ]
        if (nrow(chrom_seg) > 1) {
            chrom_seg = join_segments(chrom_seg)
        }
        chrom_seg
    })
    
    if (nrow(segs) == 0) {
        return(
            list(
                ntelomeric_ai = NA, # telomeric AI
                telomeric_ai_mean_size = NA,
                ninterstitial_ai = NA, # interstitial AI
                interstitial_ai_mean_size = NA,
                ncentromeric_ai = NA, # chromosomal AI
                ntelomeric_loh = NA, # telomeric LOH
                telomeric_loh_mean_size = NA,
                ninterstitial_loh = NA, # interstitial LOH
                interstitial_loh_mean_size = NA,
                ncentromeric_loh = NA # chromosomal LOH
            ) 
        )
    }
    
    # Add a column to call AI
    # Codes for telomeric/interstitial/whole chromosome:
    # 0 = no AI
    # 1 = telomeric
    # 2 = interstitial
    # 3 = whole chromosome
    segs$AI = NA
    segs$chrom_ploidy = NA
    
    # Loop through chromosomes
    for (chr in unique(segs$chrom)) {

        # Subset on chromosome, proceed to next if no segment
        chrom_segs = segs[which(segs$chrom == chr), ]
        if(nrow(chrom_segs) == 0) next 
        
        # Determine major ploidy for chromosome
        chrom_ploidy = group_by(chrom_segs, tcn) %>% 
            summarize(tcn_total = sum(length, na.rm = T)) %>%
            filter(tcn_total == max(tcn_total) & tcn > 0) %>% 
            pull(tcn)
        
        if (length(chrom_ploidy) == 0) { next } 
        chrom_segs$chrom_ploidy = as.integer(chrom_ploidy) # update "ploidy" column, so the new calculated value can be returned
        
        if (chrom_ploidy %% 2 == 0) { # if even
            chrom_segs$AI = c(0, 2)[match(chrom_segs$mcn == chrom_segs$lcn, c('TRUE', 'FALSE'))]
            chrom_segs$AI[which(is.na(chrom_segs$AI))] = 0 # adjust NAs
        } else if (chrom_ploidy %% 2 != 0) { # if odd
            chrom_segs$AI = c(0, 2)[match(chrom_segs$mcn + chrom_segs$lcn == ploidy &
                                              chrom_segs$lcn != 0, c('TRUE', 'FALSE'))]
            chrom_segs$AI[which(is.na(chrom_segs$AI))] = 0 # adjust NAs
        }
        
        # Put back into original seg
        segs$chrom_ploidy[which(segs$chrom == chr)] = chrom_ploidy
        segs$AI[which(segs$chrom == chr)] = chrom_segs$AI
        
        # Check relative position to centromere
        if(chrom_segs$AI[1] == 2 & nrow(chrom_segs) != 1 & chrom_segs$end[1] < (sample_chrom_info$centromere[chr])){
            segs$AI[which(segs$chrom == chr)][1] = 1 # if the first segment of chromosome is AI and does not extend to centromere --> telomeric AI
        }
        if(chrom_segs$AI[nrow(chrom_segs)] == 2 & nrow(chrom_segs) != 1 & chrom_segs$start[nrow(chrom_segs)] > (sample_chrom_info$centromere[chr])){
            segs$AI[which(segs$chrom == chr)[nrow(chrom_segs)]] = 1 # if the last segment of chromosome is AI and starts beyond the centromere --> telomeric AI
        }
        if(nrow(segs[which(segs$chrom == chr), ]) == 1 & segs$AI[which(segs$chrom == chr)][1] != 0){
            segs$AI[which(segs$chrom == chr)[1]] = 3 # if only one segment on chromosome and AI --> chromosomal AI
        }
    }
    
    # Prepare return 
    segs_loh = segs[which(segs$lcn == 0), ]
    list(
        ntelomeric_ai = nrow(segs[which(segs$AI == 1), ]), # telomeric AI
        telomeric_ai_mean_size = mean(segs$end[which(segs$AI == 1)] - segs$start[which(segs$AI == 1)]),
        ninterstitial_ai = nrow(segs[which(segs$AI == 2), ]), # interstitial AI
        interstitial_ai_mean_size = mean(segs$end[which(segs$AI == 2)] - segs$start[which(segs$AI == 2)]),
        ncentromeric_ai = nrow(segs[which(segs$AI == 3), ]), # chromosomal AI
        ntelomeric_loh = nrow(segs_loh[which(segs_loh$AI == 1), ]), # telomeric LOH
        telomeric_loh_mean_size = mean(segs_loh$end[which(segs_loh$AI == 1)] - segs_loh$start[which(segs_loh$AI == 1)]),
        ninterstitial_loh = nrow(segs_loh[which(segs_loh$AI == 2), ]), # interstitial LOH
        interstitial_loh_mean_size = mean(segs_loh$end[which(segs_loh$AI == 2)] - segs_loh$start[which(segs_loh$AI == 2)]),
        ncentromeric_loh = nrow(segs_loh[which(segs_loh$AI == 3), ]) # chromosomal LOH
    ) 
}

#' @export 
#' @rdname copy_number_scores
calculate_lst = function(segs,
                         ploidy,
                         genome = c('hg19', 'hg18', 'hg38'),
                         algorithm = c('em', 'cncf'),
                         min_size = 10e6,
                         min_probes = 50) {
    
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = FALSE)
    
    # Centromere locations
    genome = get(match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE))
    
    # Create chrom_info for sample
    segs = parse_segs(segs, algorithm) %>% 
        filter(chrom %in% 1:22)
    sample_chrom_info = get_sample_genome(segs, genome)
    
    # Count LSTs
    lst = c()
    for (chr in unique(segs$chrom)) {
        
        chrom_segs = segs[segs$chrom == chr, ]
        if (nrow(chrom_segs) < 2) next
        
        # Split into chromosome arms
        p_arm = chrom_segs[chrom_segs$start <= sample_chrom_info$centromere[sample_chrom_info$chr == chr], ] # segments starting on p arm, these might overlap centromere
        q_arm = chrom_segs[chrom_segs$end >= sample_chrom_info$centromere[sample_chrom_info$chr == chr], ] # segments ending on q arm
        p_arm = join_segments(p_arm) # shrink segments with same CN
        q_arm = join_segments(q_arm)
        if (nrow(p_arm) > 0) p_arm[nrow(p_arm), 'end'] = sample_chrom_info$centromere[sample_chrom_info$chr == chr]  # cut p-arm segment spanning centromere at centromere start
        if (nrow(q_arm) > 0) q_arm[1, 'start'] = sample_chrom_info$centromere[sample_chrom_info$chr == chr]  # set first q arm segment to start at centromere end
        
        # P arm
        # Smoothen 3-Mb segments
        n_3mb = which((p_arm$end-p_arm$start) < 3e6)
        while (length(n_3mb) > 0) {
            p_arm = p_arm[-(n_3mb[1]), ] # non-juxtaposed segments will be removed, which is fine
            p_arm = join_segments(p_arm)
            n_3mb = which((p_arm$end - p_arm$start) < 3e6)
        }
        
        # Now check for LST
        if (nrow(p_arm) >= 2) { # if more than one segment
            p_arm = cbind(p_arm, c(0, 1)[match((p_arm$end - p_arm$start) >= min_size, c('FALSE', 'TRUE'))]) # mark segments that pass length test
            for (k in 2:nrow(p_arm)) {
                if (p_arm[k, ncol(p_arm)] == 1 & p_arm[(k - 1), ncol(p_arm)] == 1 &
                    (p_arm[k, 'start'] - p_arm[(k-1), 'end']) < 3e6) { # if two juxtaposed segments are 10 Mb and the space between them is less than 3 Mb...
                    lst = c(lst, 1) # ...then add to LST
                }
            }
        }
        
        # Q arm
        # Smoothen 3-Mb segments
        n_3mb = which((q_arm$end-q_arm$start) < 3e6)
        while (length(n_3mb) > 0) {
            q_arm = q_arm[-(n_3mb[1]), ] # non-juxtaposed segments will be removed, which is fine
            q_arm = join_segments(q_arm)
            n_3mb = which((q_arm$end - q_arm$start) < 3e6)
        }
        
        # Now check for LST
        if (nrow(q_arm) >= 2) { # if more than one segment
            q_arm = cbind(q_arm, c(0, 1)[match((q_arm$end - q_arm$start) >= min_size, c('FALSE', 'TRUE'))]) # mark segments that pass length test
            for (k in 2:nrow(q_arm)) {
                if (q_arm[k, ncol(q_arm)] == 1 & q_arm[(k-1), ncol(q_arm)] == 1 &
                    (q_arm[k, 'start'] - q_arm[(k-1), 'end']) < 3e6){ # if two juxtaposed segments are 10 Mb and the space between them is less than 3 Mb...
                    lst = c(lst, 1) # ...then add to LST
                }
            }
        }
    }
    
    # Return values
    list(lst = sum(lst))
}

#' @export 
#' @rdname copy_number_scores
calculate_hrdloh = function(segs,
                            ploidy,
                            algorithm = c('em', 'cncf')) {
    
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = FALSE)
    
    segs = parse_segs(segs, algorithm) %>% 
        filter(chrom %in% 1:22)
    
    chr_del = vector()
    for (j in unique(segs$chrom)) {
        if (all(!is.na(segs$lcn[which(segs$chrom == j)]) & all(segs$lcn[which(segs$chrom == j)] == 0))) {
            chr_del = c(chr_del, j) # one parental copy of chromosome completely lost
        }
    }
    
    # Check for LOH
    segs_loh = segs[!is.na(segs$lcn), ] # remove segments with lcn = NA
    segs_loh = segs_loh[which(segs_loh$lcn == 0 & segs_loh$mcn != 0), ]
    segs_loh = segs_loh[which(segs_loh$length > 15e6), ] # check if segment is long enough
    if (length(chr_del) > 0) {
        segs_loh = segs_loh[-which(segs_loh$chrom %in% chr_del), ] # remove if whole chromosome lost
    }
    
    # Return values
    list(hrd_loh = nrow(segs_loh))
}

# Helper functions ------------------------------------------------------------------------------------------------
parse_segs = function(segs, algorithm = c('em', 'cncf')) {
    
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = FALSE)
    
    if (algorithm == 'em') {
        segs$tcn = segs$tcn.em
        segs$lcn = segs$lcn.em
        segs$cf = segs$cf.em
    }
    mutate(segs,
           length = end - start,
           lcn = ifelse(tcn <= 1, 0, lcn),  # correct mcn, lcn for cases of tcn = 1 // sometimes FACETS set lcn = NA when tcn = 1, when it clearly has to be 0
           mcn = tcn - lcn) %>% 
        select(-ends_with('\\.em'))
}

get_sample_genome = function(segs, genome) {
    map_dfr(unique(segs$chrom), function(x) {
        centromere = genome$centromere[which(genome$chrom == x)]
        chrstart = min(segs$start[which(segs$chrom == x)])
        chrend = max(segs$end[which(segs$chrom == x)])
        size = chrend - chrstart
        plength = centromere - chrstart
        if (plength < 0) plength = 0 # acrocentric chromosomes
        qlength = chrend - centromere
        data.frame(
            chr = x,
            centromere = centromere,
            chrstart = chrstart,
            chrend = chrend,
            size = size,
            plength = plength,
            qlength = qlength
        )
    }
    )
}

is_genome_doubled = function(segs, chrom_info, treshold = 0.5) {
    
    # Calculated length of autosomal
    autosomal_genome = sum(as.numeric(chrom_info$size[chrom_info$chr %in% 1:22]))
    
    # Check for whole-genome duplication // PMID 30013179
    frac_elevated_mcn = sum(as.numeric(segs$length[which(segs$mcn >= 2 & segs$chrom %in% 1:22)])) / autosomal_genome
    wgd = frac_elevated_mcn > treshold
    
    wgd
}

# Join segments with identical copy number that are separated for some reason
# Logic: if different tcn ==> don't join
# if any NAs in mcn/lcn ==> don't join
# if identical mcn/lcn content ==> do join
join_segments = function(chrom_seg) {
    if (nrow(chrom_seg) < 2) {
        chrom_seg
    } else {
        new_chr = chrom_seg
        seg_class = c(1)
        for (j in 2:nrow(new_chr)) {
            # if adjacent segments have same allelic content, assign to same class
            if (new_chr[(j - 1), 'tcn'] != new_chr[j, 'tcn']) { # if tcn differs, definitely don't condense
                seg_class = c(seg_class, seg_class[j - 1] + 1)
            } else if (any(is.na(c(new_chr[(j - 1), 'mcn'], new_chr[(j), 'mcn'])))) { # but if tcn the same, check for difference in allelic content
                seg_class = c(seg_class, seg_class[j - 1] + 1)
            } else if (new_chr[(j - 1), 'mcn'] == new_chr[j, 'mcn'] &
                       new_chr[(j - 1), 'lcn'] == new_chr[j, 'lcn']) {
                seg_class = c(seg_class, seg_class[j - 1])
            } else {
                seg_class = c(seg_class, seg_class[j - 1] + 1)
            }
        }
        for (j in unique(seg_class)) {
            # condense segments belonging to same class
            new_chr[seg_class %in% j, 'end'] = max(new_chr[seg_class %in% j, 'end'])
            new_chr[seg_class %in% j, 'num.mark'] = sum(new_chr[seg_class %in% j, 'num.mark'])
        }
        new_chr[!duplicated(seg_class), ]
    }
}
