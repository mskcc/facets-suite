#' Copy-number based scores
#' 
#' Calculate the following 
#' \itemize{
#'  \item{Copy-number statistics:} {Fraction of genome altered, altered chromosome arms and genome doubling flag.}
#'  \item{LST score:} {Large-scale state transitions, see source URL.}
#'  \item{NtAI:} {Telomeric allelic imbalance, see source URL.}
#'  \item{HRD-LOH:} {HRD-LOH score, see source URL.}
#' }
#'
#' @param segs FACETS segmentation output.
#' @param ploidy Sample ploidy.
#' @param genome Genome build.
#' @param algorithm Choice between FACETS \code{em} and \code{cncf} algorithm.
#' @param min_size Minimum length of segment, as defined per function.
#' @param min_probes Minimum number of SNPs per segment, as defined per function.
#' 
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/29622463}
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/26015868}
#'
#' @return `data.frame` with one or more values from function.
#'
#' @import dplyr
#' @importFrom purrr map_dfr map_lgl discard
#' 
#' @name copy_number_scores
NULL

#' @export 
#' @rdname copy_number_scores
calculate_fcna = function(segs,
                          ploidy,
                          genome = c('hg19', 'hg18', 'hg38'),
                          algorithm = c('em', 'cncf')) {
    
    genome = match.arg(genome)
    algorithm = match.arg(algorithm)

    # Centromere locations
    genome = get(genome)

    # Create chrom_info for sample
    segs = parse_segs(segs, algorithm)
    sample_chrom_info = get_sample_genome(segs, genome)

    # Calculated length of interrogated genome
    interrogated_genome = sum(sample_chrom_info$size)
    autosomal_genome = sum(sample_chrom_info$size[sample_chrom_info$chr %in% 1:22])

    # Check for whole-genome duplication // PMID 30013179
    wgd_treshold = 0.5 # treshold
    frac_elevated_mcn = sum(segs$length[which(segs$mcn >= 2 & segs$chrom %in% 1:22)]) / autosomal_genome
    wgd = frac_elevated_mcn > wgd_treshold

    # Calculate fraction of genome altered
    sample_ploidy = ifelse(wgd, round(ploidy), 2)
    if (!wgd) {
        diploid_length = sum(segs$length[which(segs$tcn == sample_ploidy & segs$lcn == 1)])
    } else if (wgd) {
        diploid_length = sum(segs$length[which(segs$tcn == sample_ploidy & segs$lcn >= 1)])
    }
    frac_altered = (interrogated_genome - diploid_length) / interrogated_genome

    # Altered arms, slightly adjusted to PMID 29622463
    segs = left_join(segs, select(sample_chrom_info, chr, centromere), by = c('chrom' = 'chr'))
    altered_arms = sapply(unique(segs$chrom), function(x) {
        seg_p = segs[which(segs$chrom == x & segs$start < segs$centromere), ]
        seg_p$end = ifelse(seg_p$end > seg_p$centromere, seg_p$centromere, seg_p$end)
        seg_p$length = seg_p$end - seg_p$start
        seg_q = segs[which(segs$chrom == x & segs$end > segs$centromere), ]
        seg_q$start = ifelse(seg_q$start < seg_q$centromere, seg_q$centromere, seg_q$start)
        seg_q$length = seg_q$end - seg_q$start

        if (!wgd) {
            seg_p_unaltered = sum(seg_p$length[which(seg_p$tcn == sample_ploidy & seg_p$lcn == 1)])
            seg_q_unaltered = sum(seg_q$length[which(seg_q$tcn == sample_ploidy & seg_q$lcn == 1)])
        } else if (wgd) {
            seg_p_unaltered = sum(seg_p$length[which(seg_p$tcn == sample_ploidy & seg_p$lcn >= 1)])
            seg_q_unaltered = sum(seg_q$length[which(seg_q$tcn == sample_ploidy & seg_q$lcn >= 1)])
        }
        paste0(paste0(x, c('p', 'q'))[c((sample_chrom_info$plength[sample_chrom_info$chr == x] - 
                                             seg_p_unaltered) / sample_chrom_info$plength[sample_chrom_info$chr == x] > .8,
                                        (sample_chrom_info$qlength[sample_chrom_info$chr == x] -
                                             seg_q_unaltered) / sample_chrom_info$qlength[sample_chrom_info$chr == x] > .8) %>%
                                          map_lgl(~ifelse(is.na(.), FALSE, .))], collapse = ',')
    }
    ) %>% strsplit(., ',') %>%
        unlist %>%
        discard(. %in% c('', '13p', '14p', '15p', '21p', '22p'))

    # weighted fraction copy-number altered
    frac_altered_w = select(sample_chrom_info, chr, p = plength, q = qlength) %>%
        gather(arm, length, -chr) %>%
        filter(!paste0(chr, arm) %in% c('13p', '14p', '15p', '21p', '22p')) %>%
        summarize(sum(length[paste0(chr, arm) %in% altered_arms]) / sum(length)) %>%
        as.numeric()

    data.frame(
        genome_doubled = wgd,
        fcna = frac_altered,
        weighted_fcna = frac_altered_w,
        aneuploidy_score = length(altered_arms),
        altered_arms = paste0(altered_arms, collapse = ',')
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
    
    genome = match.arg(genome)
    algorithm = match.arg(algorithm)
    
    # Centromere locations
    genome = get(genome)
    
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
        joined_segs = rbind(joined_segs, chrom_seg)
    })
    
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
            summarize(tcn_total = sum(length)) %>%
            filter(tcn_total == max(tcn_total) & tcn > 0) %>% 
            pull(tcn)
        chrom_segs$chrom_ploidy = chrom_ploidy # update "ploidy" column, so the new calculated value can be returned
        
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
    data.frame(
        ntelomeric_ai = nrow(segs[which(segs$AI == 1), ]), # telomeric AI
        telomeric_ai_mean_size = mean(segs$end[which(segs$AI == 1)] - segs$start[which(segs$AI == 1)]),
        ninterstitial_ai = nrow(segs[which(segs$AI == 2), ]), # interstitial AI
        interstitial_ai_mean_size = mean(segs$end[which(segs$AI == 2)] - segs$start[which(segs$AI == 2)]),
        ncentromeric_ai = nrow(segs[which(segs$AI == 3), ]), # chromosomal AI
        ntelomeric_loh = nrow(segs_loh[which(segs_loh$AI == 1), ]), # telomeric LOH
        telomeric_loh_mean_size = mean(segs_loh$end[which(segs_loh$AI == 1)] - segs_loh$start[which(segs_loh$AI == 1)]),
        ninterstitial_loh = nrow(seg_loh[which(seg_loh$AI == 2), ]), # interstitial LOH
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
    
    genome = match.arg(genome)
    algorithm = match.arg(algorithm)
    
    # Centromere locations
    genome = get(genome)
    
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
        while(length(n_3mb) > 0) {
            p_arm = p_arm[-(n_3mb[1]), ] # non-juxtaposed segments will be removed, which is fine
            p_arm = join_segments(p_arm)
            n_3mb = which((p_arm$end - p_arm$start) < 3e6)
        }
        
        # Now check for LST
        if (nrow(p_arm) >= 2) { # if more than one segment
            p_arm = cbind(p_arm, c(0, 1)[match((p_arm$end - p_arm$start) >= min_length, c('FALSE', 'TRUE'))]) # mark segments that pass length test
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
        while(length(n_3mb) > 0) {
            q_arm = q_arm[-(n_3mb[1]), ] # non-juxtaposed segments will be removed, which is fine
            q_arm = join_segments(q_arm)
            n_3mb = which((q_arm$end - q_arm$start) < 3e6)
        }
        
        # Now check for LST
        if (nrow(q_arm) >= 2) { # if more than one segment
            q_arm = cbind(q_arm, c(0, 1)[match((q_arm$end - q_arm$start) >= min_length, c('FALSE', 'TRUE'))]) # mark segments that pass length test
            for (k in 2:nrow(q_arm)) {
                if (q_arm[k, ncol(q_arm)] == 1 & q_arm[(k-1), ncol(q_arm)] == 1 &
                    (q_arm[k, 'start'] - q_arm[(k-1), 'end']) < 3e6){ # if two juxtaposed segments are 10 Mb and the space between them is less than 3 Mb...
                    lst = c(lst, 1) # ...then add to LST
                }
            }
        }
    }
    
    # Return values
    data.frame(lst = lst)
}

#' @export 
#' @rdname copy_number_scores
calculate_hrdloh = function(segs,
                            ploidy,
                            genome = c('hg19', 'hg18', 'hg38'),
                            algorithm = c('em', 'cncf')) {
    
    segs = filter(segs, chrom %in% 1:22)
    
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
    segs_loh = segs_loh[!which(segs_loh$chrom %in% chr_del), ] # remove if whole chromosome lost
    
    # Return values
    data.frame(hrd_loh = nrow(segs_loh))
}

# Helper functions ------------------------------------------------------------------------------------------------
parse_segs = function(segs, algorithm = c('em', 'cncf')) {
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
        size = chrend - chrstart
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
        for(j in 2:nrow(new_chr)) {
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
        for(j in unique(seg_class)) {
            # condense segments belonging to same class
            new_chr[seg_class %in% j, 'end'] = max(new_chr[seg_class %in% j, 'end'])
            new_chr[seg_class %in% j, 'num.mark'] = sum(new_chr[seg_class %in% j, 'num.mark'])
        }
        new_chr[!duplicated(seg_class), ]
    }
}
