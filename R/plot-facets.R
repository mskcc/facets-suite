#' Plot FACETS output
#' 
#' Generate the following plots:
#' \itemize{
#'  \item{Log-ratio:} {Copy-number segmentation based on tumor-normal read coverage comparison.}
#'  \item{Allelic imbalance:} {Allelic imbalance based on somatic changes in zygosity of heterozygous SNPs in the normal.}
#'  \item{Integer copy number:} {Inference of the major and minor copy-number states based on the tumor-normal log ratio and allelic content.}
#'  \item{Cellular fraction:} {Estimate of fraction of cells in sample harboring each segment.}
#' }
#'
#' @param facets_data Output object from \code{run_facets}.
#' @param colors Vector of two colors for alternating chromosomes.
#' @param plotX If \code{TRUE}, includes chromosome X.
#' @param genome Genome build.
#' @param highlight_gene Highlight gene(s), provide gene symbol or mapped positions (internally).
#' @param adjust_diplogr Normalize by sample dipLogR.
#' @param method When available, choose between plotting solution from \code{em} or \code{cncf} algorithm.
#' @param subset_snps Subset the SNP profile to reduce weight of plotting, supply a factor by which to reduce or \code{TRUE} for default.
#' @param plot_chroms Chromosomes to plot when using \code{closeup_plot}.
#' @param return_object If \code{TRUE}, returns \code{ggplot2} object instead of printing plot.
#' 
#' @return \code{ggplot2} objects, see input parameter \code{return_object}.
#' 
#' @import ggplot2
#' @importFrom dplyr filter mutate
#' @importFrom grDevices colorRampPalette
#' @importFrom egg ggarrange
#' @importFrom scales pretty_breaks
#'
#' @name plot_facets
NULL

#' @export
#' @rdname plot_facets
cnlr_plot = function(facets_data,
                     colors = c('#0080FF', '#4CC4FF'),
                     plotX = FALSE,
                     genome = c('hg19', 'hg18', 'hg38'),
                     highlight_gene = NULL,
                     adjust_diplogr = TRUE,
                     subset_snps = NULL,
                     return_object = FALSE) {
    
    genome = match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE)
    
    snps = facets_data$snps
    segs = facets_data$segs
    diplogr = facets_data$diplogr
    if (!plotX) {
        snps = subset(snps, chrom < 23)
        segs = subset(segs, chrom < 23)
    }
    
    snps = get_cum_chr_maploc(snps, genome)
    mid = snps$mid[names(snps$mid) %in% snps$snps$chrom]
    centromeres = snps$centromeres
    snps = snps$snps
    
    snps$cnlr_median = rep(segs$cnlr.median, segs$num.mark)
    
    starts = cumsum(c(1, segs$num.mark))[seq_along(segs$num.mark)]
    ends = cumsum(c(segs$num.mark))
    my_starts = snps[starts, c('chr_maploc', 'cnlr_median')]
    my_ends = snps[ends, c('chr_maploc', 'cnlr_median')]
    
    ymin = floor(min(segs$cnlr.median, na.rm = T))
    ymax = ceiling(max(segs$cnlr.median, na.rm = T))
    if (ymin > -3) ymin = -3
    if (ymax < 3) ymax = 3
    
    if (adjust_diplogr) {
        snps$cnlr = snps$cnlr - diplogr
        my_starts$cnlr_median = my_starts$cnlr_median - diplogr
        my_ends$cnlr_median = my_ends$cnlr_median - diplogr
        diplogr = Inf
    }
    
    if (!is.null(subset_snps)) {
        if (subset_snps == TRUE) {
            snps = subset_snps(snps)
        } else if (is.numeric(subset_snps)) {
            snps = subset_snps(snps, subset_snps)
        }
    }
    
    pt_cols = colors[c(snps$chrom %% 2) + 1]
    
    # plot
    cnlr = ggplot(snps) +
        geom_point(aes(y = cnlr, x = chr_maploc), pch = 19, col = pt_cols, size = .4) +
        scale_x_continuous(breaks = mid, labels = names(mid), expand = c(.01, 0)) +
        scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(ymin, ymax)) +
        geom_hline(yintercept = diplogr, color = 'sandybrown', size = .8) +
        geom_segment(data = segs, aes(x = my_starts$chr_maploc, xend = my_ends$chr_maploc,
                                      y = my_starts$cnlr_median, yend = my_ends$cnlr_median),
                     col = 'red3', size = 1, lineend = 'butt') +
        labs(x = NULL, y = 'Copy number log ratio') +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 0, size = 8, color = 'black'),
              axis.text.y = element_text(angle = 0, size = 8, color = 'black'),
              text = element_text(size = 10),
              panel.grid.minor.x = element_line(colour = 'grey', size = .2),
              panel.grid.major.x = element_line(colour = 'grey', size = 0),
              plot.margin = unit(c(0, 1, 0, 0), 'lines'))
    
    # if highligthing gene
    if (!is.null(highlight_gene)) {
        if (is.character(highlight_gene)) {
            highlight_gene = get_gene_position(highlight_gene)
        }
        snps$gene = FALSE
        snps$gene[which(snps$chrom %in% highlight_gene$chrom &
                        snps$chr_maploc >= highlight_gene$start &
                        snps$chr_maploc <= highlight_gene$end)] = TRUE
        cnlr = cnlr +
            geom_vline(xintercept = highlight_gene$mid, color = 'palevioletred1') +
            geom_point(data = filter(snps, gene == TRUE), aes(y = cnlr, x = chr_maploc), color = '#525252', size = .4) 
    }
    
    if (return_object == TRUE) {
        cnlr 
    } else {
        suppressMessages(print(cnlr))
    }
}

#' @export
#' @rdname plot_facets
valor_plot = function(facets_data,
                      colors = c('#0080FF', '#4CC4FF'),
                      plotX = FALSE,
                      genome = c('hg19', 'hg18', 'hg38'),
                      highlight_gene = NULL,
                      subset_snps = NULL,
                      return_object = FALSE) {
    
    genome = match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE)
    
    snps = facets_data$snps
    segs = facets_data$segs
    diplogr = facets_data$diplogr
    if (!plotX) {
        snps = subset(snps, chrom < 23)
        segs = subset(segs, chrom < 23)
    }
    
    snps = get_cum_chr_maploc(snps, genome)
    mid = snps$mid[names(snps$mid) %in% snps$snps$chrom]
    centromeres = snps$centromeres
    snps = snps$snps
    
    snps$mafr = rep(sqrt(abs(segs$mafR)), segs$num.mark)
    
    starts = cumsum(c(1, segs$num.mark))[seq_along(segs$num.mark)]
    ends = cumsum(c(segs$num.mark))
    my_starts = snps[starts, c('chr_maploc', 'mafr')]
    my_ends = snps[ends, c('chr_maploc', 'mafr')]
    
    if (!is.null(subset_snps)) {
        if (subset_snps == TRUE) {
            snps = subset_snps(snps)
        } else if (is.numeric(subset_snps)) {
            snps = subset_snps(snps, subset_snps)
        }
    }
    
    pt_cols = colors[c(snps$chrom %% 2) + 1]
    
    # plot
    valor = ggplot(snps) +
        geom_point(aes(y = valor, x = chr_maploc), pch = 19, col = pt_cols, size = .4) +
        scale_x_continuous(breaks = mid, labels = names(mid), expand = c(.01, 0)) +
        geom_segment(data = segs, col = 'red3', size = 1, 
                     aes(x = my_starts$chr_maploc, xend = my_ends$chr_maploc,
                         y = my_starts$mafr, yend = my_ends$mafr)) +
        geom_segment(data = segs, col = 'red3', size = 1, 
                     aes(x = my_starts$chr_maploc, xend = my_ends$chr_maploc,
                         y = -my_starts$mafr, yend = -my_ends$mafr)) +
        labs(x = NULL, y = 'Variant allele log odds ratio') +
        ylim(-4, 4) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 0, size = 8, color = 'black'),
              axis.text.y = element_text(angle = 0, size = 8, color = 'black'),
              text = element_text(size = 10),
              panel.grid.minor.x = element_line(colour = 'grey', size = .2),
              panel.grid.major.x = element_line(colour = 'grey', size = 0),
              plot.margin = unit(c(0, 1, 0, 0), 'lines'))
    
    # if highligthing gene, can take gene name(s) or already mapped positions
    if (!is.null(highlight_gene)) {
        if (is.character(highlight_gene)) {
            gene_pos = get_gene_position(gene_pos)
        } 
        snps$gene = FALSE
        snps$gene[which(snps$chrom %in% highlight_gene$chrom &
                        snps$chr_maploc >= highlight_gene$start &
                        snps$chr_maploc <= highlight_gene$end)] = TRUE
        valor = valor +
            geom_vline(xintercept = highlight_gene$mid, color = 'palevioletred1') +
            geom_point(data = filter(snps, gene == TRUE), aes(y = cnlr, x = chr_maploc), color = '#525252', size = .4) 
    }
    
    if (return_object == TRUE) {
        valor 
    } else {
        suppressMessages(print(valor))
    }
}

#' @export
#' @rdname plot_facets
cf_plot = function(facets_data,
                   method = c('em', 'cncf'),
                   plotX = FALSE,
                   genome = c('hg19', 'hg18', 'hg38'),
                   return_object = FALSE) {
    
    genome = match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE)
    
    snps = facets_data$snps
    segs = facets_data$segs
    if (!plotX) {
        snps = subset(snps, chrom < 23)
        segs = subset(segs, chrom < 23)
    }
    
    snps = get_cum_chr_maploc(snps, genome)
    mid = snps$mid[names(snps$mid) %in% snps$snps$chrom]
    centromeres = snps$centromeres
    snps = snps$snps
    
    if (method == 'em') {
        cols = c(grDevices::colorRampPalette(c('white', 'steelblue'))(10), 'papayawhip')[round(10 * segs$cf.em + 0.501)]
        my_ylab = 'CF (EM)'
    } else if (method == 'cncf') {
        my_ylab = 'CF (CNCF)'
        cols = c(grDevices::colorRampPalette(c('white', 'steelblue'))(10), 'papayawhip')[round(10 * segs$cf + 0.501)]
    }
    
    starts = cumsum(c(1, segs$num.mark))[seq_along(segs$num.mark)]
    ends = cumsum(c(segs$num.mark))
    
    cf = ggplot(segs) +
        geom_rect(aes(xmin = starts, xmax = ends, ymax = 1, ymin = 0),
                  fill = cols, col = 'white', size = 0) +
        scale_x_continuous(breaks = mid, labels = names(mid), expand = c(.01, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(x = NULL, y = my_ylab) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 0, size = 8, color = 'black'),
              axis.text.y = element_text(angle = 0, size = 8, color = 'white'),
              axis.ticks.y = element_line(color = 'white'),
              text = element_text(size = 10),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              plot.margin = unit(c(0, 1, 0, 0), 'lines'))
    
    if (return_object == TRUE) {
        cf 
    } else {
        suppressMessages(print(cf))
    }
}

#' @export
#' @rdname plot_facets
icn_plot = function(facets_data,
                    method = c('em', 'cncf'),
                    plotX = FALSE,
                    genome = c('hg19', 'hg18', 'hg38'),
                    highlight_gene = NULL,
                    return_object = FALSE) {
    
    genome = match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE)
    
    snps = facets_data$snps
    segs = facets_data$segs
    if (!plotX) {
        snps = subset(snps, chrom < 23)
        segs = subset(segs, chrom < 23)
    }
    
    snps = get_cum_chr_maploc(snps, genome)
    mid = snps$mid[names(snps$mid) %in% snps$snps$chrom]
    centromeres = snps$centromeres
    snps = snps$snps
    
    if (method == 'em') {
        tcnscaled = segs$tcn.em
        tcnscaled[segs$tcn.em > 5 & !is.na(segs$tcn.em)] = (5 + (tcnscaled[segs$tcn.em > 5 & !is.na(segs$tcn.em)] - 5) / 3)
        lcn = rep(segs$lcn.em, segs$num.mark)
        my_ylab = 'Integer copy number (EM)'        
    } else if (method == 'cncf') {
        tcnscaled = segs$tcn
        tcnscaled[segs$tcn > 5 & !is.na(segs$tcn)] = (5 + (tcnscaled[segs$tcn > 5 & !is.na(segs$tcn)] - 5) / 3)
        lcn = rep(segs$lcn, segs$num.mark)
        my_ylab = 'Integer copy number (CNCF)'        
    }
    tcn = rep(tcnscaled, segs$num.mark)

    snps$tcn = tcn
    snps$lcn = lcn
    starts = cumsum(c(1, segs$num.mark))[seq_along(segs$num.mark)]
    ends = cumsum(c(segs$num.mark))
    my_tcn_starts = snps[starts, c('chr_maploc', 'tcn')]
    my_tcn_ends = snps[ends, c('chr_maploc', 'tcn')]
    my_lcn_starts = snps[starts, c('chr_maploc', 'lcn')]
    my_lcn_ends = snps[ends, c('chr_maploc', 'lcn')]
    
    icn = ggplot(segs) +
        geom_segment(col = 'red', size = 1, 
                     aes(x = my_lcn_starts$chr_maploc, xend = my_lcn_ends$chr_maploc, 
                         y = my_lcn_starts$lcn, yend = my_lcn_ends$lcn)) +
        geom_segment(col = 'black', size = 1, 
                     aes(x = my_tcn_starts$chr_maploc, xend = my_tcn_ends$chr_maploc, 
                         y = my_tcn_starts$tcn, yend = my_tcn_ends$tcn)) +
        scale_y_continuous(breaks=c(0:5, 5 + seq_len(35) / 3), labels = 0:40, limits = c(0, NA)) +
        scale_x_continuous(breaks = mid, labels = names(mid), expand = c(.01, 0)) +
        labs(x = NULL, y = my_ylab) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 0, size = 8, color = 'black'),
              axis.text.y = element_text(angle = 0, size = 8, color = 'black'),
              text = element_text(size = 10),
              panel.grid.minor.x = element_line(colour = 'grey', size = .2),
              panel.grid.major.x = element_line(colour = 'grey', size = 0),
              plot.margin = unit(c(0, 1, 0, 0), 'lines'))

    if (!is.null(highlight_gene)) {
        if (is.character(highlight_gene)) {
            gene_pos = get_gene_position(gene_pos)
        }
        icn = icn +
            geom_vline(xintercept = gene_pos$mid, color = 'palevioletred1')
    }
    
    if (return_object == TRUE) {
        icn 
    } else {
        suppressMessages(print(icn))
    }
}

#' @export
#' @rdname plot_facets
closeup_plot = function(facets_data,
                        highlight_gene = NULL,
                        plot_chroms = NULL,
                        method = c('em', 'cncf'),
                        genome = c('hg19', 'hg18', 'hg38'),
                        return_object = FALSE) {
    
    genome = match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE)
    method = match.arg(method, c('em', 'cncf'), several.ok = FALSE)
    
    # If highlighting a gene
    if (!is.null(highlight_gene)) {
        highlight_gene = get_gene_position(highlight_gene, genome)
    }
    
    # If no chromosome is chose, gene has to be 
    if (is.null(plot_chroms)) {
        if (is.null(highlight_gene)) stop('Specify at least a gene or chromosome to zoom in on.', call. = F)
        if (nrow(highlight_gene) > 1) {
            plot_chroms = min(highlight_gene$chrom):max(highlight_gene$chrom)
        } else {
            plot_chroms = (highlight_gene$chrom-1):(highlight_gene$chrom+1)
        }
    }
    
    facets_data$snps = filter(facets_data$snps, chrom %in% plot_chroms)
    facets_data$segs = filter(facets_data$segs, chrom %in% plot_chroms)

    cnlr = cnlr_plot(facets_data, genome = genome, highlight_gene = highlight_gene, return_object = TRUE)
    valor = valor_plot(facets_data, genome = genome, highlight_gene = highlight_gene, return_object = TRUE)
    icn = icn_plot(facets_data, genome = genome, method = method, highlight_gene = highlight_gene, return_object = T)
    
    if (return_oject == TRUE) {
        list(cnlr, valor, icn)
    } else {
        suppressWarnings(ggarrange(plots = list(cnlr, valor, icn), ncol = 1))
    }
    
}

# Helper functions ------------------------------------------------------------------------------------------------
# Get positions of snps, running accros whole genome
get_cum_chr_maploc = function(snps,
                              genome = c('hg19', 'hg18', 'hg38')) {
    
    genome = get(genome)
    
    cum_chrom_lengths = cumsum(as.numeric(genome$size))
    mid = cum_chrom_lengths - (genome$size / 2)
    names(mid) = seq_len(nrow(genome))
    centromeres = genome$centromere + c(0, cum_chrom_lengths[-length(cum_chrom_lengths)])
    
    snps$chr_maploc = snps$maploc + c(0, cum_chrom_lengths)[snps$chrom]
    
    list(snps = snps, mid = mid, centromeres = centromeres)
}

# Look up gene position in internal data tables
get_gene_position = function(genes,
                             genome = c('hg19', 'hg18', 'hg38')) {
    
    gene_loci = get(paste0('genes_', genome)) 
    genome = get(genome)
    cum_chrom_lengths = cumsum(as.numeric(genome$size))
    
    filter(gene_loci, gene %in% genes) %>%
        mutate(chrom = as.numeric(chrom),
               mid = start + (end-start) / 2,
               start = start + cum_chrom_lengths[chrom-1],
               end = end + cum_chrom_lengths[chrom-1],
               mid = mid + cum_chrom_lengths[chrom-1]) 
}

# Proproptionally subset SPNs per chromosome to a given fraction
subset_snps = function(snps, by_factor = 5) {
    set.seed(42)
    group_by(snps, chrom) %>% 
        sample_frac(1 / by_factor, replace = FALSE)
}
