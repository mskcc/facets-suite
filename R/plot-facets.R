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
#' @param cols Vector of two colors for alternating chromosomes.
#' @param subset_indices ???.
#' @param plotX Plot chromosome X.
#' @param genome Genome build.
#' @param gene_pos Highlight gene.
#' @param adjust_diplogr Normalize by sample dipLogR.
#' @param method When available, choose between plotting solution from \code{em} or \code{cncf} algorithm.
#' 
#' @return \code{ggplot2} objects.
#' 
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#'
#' @name plot_facets
NULL

#' @export
#' @rdname plot_facets
cnlr_plot = function(facets_data,
                     cols = c('#0080FF', '#4CC4FF'),
                     subset_indices = NULL,
                     plotX = FALSE,
                     genome = c('hg19', 'hg18', 'hg38'),
                     gene_pos = NULL,
                     adjust_diplogr = TRUE) {
    
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
    
    starts = cumsum(c(1, segs$num.mark))[1:length(segs$num.mark)]
    ends = cumsum(c(segs$num.mark))
    my_starts = snps[starts, c('chr_maploc', 'cnlr_median')]
    my_ends = snps[ends, c('chr_maploc', 'cnlr_median')]
    
    if (!is.null(subset_indices)) {
        snps = mat[subset_indices, ]
    }
    
    pt_cols = cols[c(snps$chrom %% 2) + 1]
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
    if (!is.null(gene_pos)) {
        snps$gene = FALSE
        snps$gene[which(snps$chrom == gene.pos$chrom &
                            snps$maploc >= gene_pos$start &
                            snps$maploc <= gene_pos$end)] = TRUE
        cnlr = cnlr + geom_vline(xintercept = gene_pos$mid, color = 'palevioletred1') +
            geom_point(data = subset(snps, gene == T), aes(y = cnlr, x = chr_maploc), color = '#525252', size = .4) 
    } else {
        cnlr   
    }
}

#' @export
#' @rdname plot_facets
valor_plot = function(facets_data,
                      cols = c('#0080FF', '#4CC4FF'),
                      subset_indices = NULL,
                      plotX = FALSE,
                      genome = c('hg19', 'hg18', 'hg38'),
                      gene_pos = NULL) {
    
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
    
    starts = cumsum(c(1, segs$num.mark))[1:length(segs$num.mark)]
    ends = cumsum(c(segs$num.mark))
    my_starts = snps[starts, c('chr_maploc', 'mafr')]
    my_ends = snps[ends, c('chr_maploc', 'mafr')]
    
    if (!is.null(subset_indices)) {
        snps = mat[subset_indices, ]
    }
    
    pt_cols = cols[c(snps$chrom %% 2) + 1]
    
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
    
    # if highligthing gene
    if (!is.null(gene_pos)) {
        snps$gene = FALSE
        snps$gene[which(snps$chrom == gene.pos$chrom &
                            snps$maploc >= gene_pos$start &
                            snps$maploc <= gene_pos$end)] = TRUE
        valor + geom_vline(xintercept = gene_pos$mid, color = 'palevioletred1') +
            geom_point(data = subset(snps, gene == T), aes(y = cnlr, x = chr_maploc), color = '#525252', size = .4) 
    } else {
        valor   
    }
}

#' @export
#' @rdname plot_facets
cf_plot = function(facets_data,
                   method = c('em', 'cncf'),
                   plotX = FALSE,
                   genome = c('hg19', 'hg18', 'hg38')) {
    
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
    
    starts = cumsum(c(1, segs$num.mark))[1:length(segs$num.mark)]
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
    cf
}

#' @export
#' @rdname plot_facets
icn_plot = function(facets_data,
                    method = c('em', 'cncf'),
                    plotX = FALSE,
                    genome = c('hg19', 'hg18', 'hg38'),
                    gene_pos = NULL) {
    
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
    
    if(method == 'em') {
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
    starts = cumsum(c(1, segs$num.mark))[1:length(segs$num.mark)]
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
        scale_y_continuous(breaks=c(0:5, 5 + (1:35) / 3), labels = 0:40, limits = c(0, NA)) +
        scale_x_continuous(breaks = mid, labels = names(mid), expand = c(.01, 0)) +
        labs(x = NULL, y = my_ylab) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 0, size = 8, color = 'black'),
              axis.text.y = element_text(angle = 0, size = 8, color = 'black'),
              text = element_text(size = 10),
              panel.grid.minor.x = element_line(colour = 'grey', size = .2),
              panel.grid.major.x = element_line(colour = 'grey', size = 0),
              plot.margin = unit(c(0, 1, 0, 0), 'lines'))

    if(!is.null(gene_pos)) {
        icn + geom_vline(xintercept = gene_pos$mid, color = 'palevioletred1')
    } else {
        icn    
    }
}

get_cum_chr_maploc = function(snps, genome = c('hg19', 'hg18', 'hg38')) {
    
    genome = get(genome)
    
    cum_chrom_lengths = cumsum(as.numeric(genome$size))
    mid = cum_chrom_lengths - (genome$size / 2)
    names(mid) = 1:nrow(genome)
    centromeres = genome$centromere + c(0, cum_chrom_lengths[-length(cum_chrom_lengths)])
    
    snps$chr_maploc = snps$maploc + c(0, cum_chrom_lengths)[snps$chrom]
    
    list(snps = snps, mid = mid, centromeres = centromeres)
}
