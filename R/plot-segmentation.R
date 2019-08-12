#' Plot segmentation profile
#' 
#' Creates an IGV-like graphical representation of the copy-number segments across the samples in a segmentation object, as output by \code{run_facets}.
#'
#' @param segs FACETS segment output, IGV formatted. 
#' @param plotX If \code{TRUE}, includes chromosome X.
#' @param sample_order Manual order of samples.
#' @param cap_log_ratios Cap log-ratios at the absolute value.
#' @param colors Vector of three colors, giving the low-, mid- and high-point of the color scale.
#' @param return_object If \code{TRUE}, returns \code{ggplot2} object instead of printing plot.
#'
#' @return \code{ggplot2} objects, see input parameter \code{return_object}.
#'
#' @importFrom dplyr mutate left_join
#' @importFrom purrr map_dfr
#' @import ggplot2
 
#' @export
plot_segmentation = function(segs,
                             plotX = FALSE,
                             sample_order = NULL,
                             cap_log_ratios = TRUE,
                             colors = c('darkblue', 'white', 'red'),
                             return_object = FALSE) {
    
    if (!plotX) { 
        segs = filter(segs, chrom != 23)
        chrom_order = seq(1:23)
    } else {
        chrom_order = seq(1:22)
    }
    segs = mutate(segs, chrom = factor(chrom, levels = chrom_order, ordered = T))
    
    if (!is.null(sample_order)) {
        if (!all(segs$ID %in% sample_order)) stop('Samples missing from provided sample order', call. = FALSE)
        segs = mutate(segs, ID = factor(ID, levels = sample_order, ordered = T))
    }
    
    # Determine lengths of chromosomes and adjust X-axis accordingly
    chrom_lengths = map_dfr(unique(segs$chrom), function(x) {
        chrom_max = max(segs$loc.end[segs$chrom == x], na.rm = T)
        chrom_min = min(segs$loc.end[segs$chrom == x], na.rm = T)
        list(chrom = x,
             chrom_length = as.numeric(chrom_max - chrom_min))
    }) %>% 
        mutate(.data, rel_length = chrom_length/sum(chrom_length))
    segs = left_join(segs, chrom_lengths, by = 'chrom')
    
    
    # Cap log-ratios and set colorlegend
    if (cap_log_ratios != FALSE) {
        if (is.numeric(cap_log_ratios)) {
            segs$seg.mean[which(segs$seg.mean > cap_log_ratios)] = cap_log_ratios
            segs$seg.mean[which(segs$seg.mean < -cap_log_ratios)] = -cap_log_ratios
            legend_breaks = c(-cap_log_ratios, -cap_log_ratios/2, 0, cap_log_ratios/2, cap_log_ratios)
        } else {
            segs$seg.mean[which(segs$seg.mean > 2)] = 2
            segs$seg.mean[which(segs$seg.mean < -2)] = -2    
            legend_breaks = seq(-2, 2, 1)
        }
    } else {
        legend_breaks = c(min(segs$seg.mean), min(segs$seg.mean)/2, 0, max(segs$seg.mean)/2, max(segs$seg.mean))
    }
    
    # Set Y-axis
    sample_number = length(unique(segs$ID))
    increments = 100 / sample_number
    max_vec = cumsum(rep(increments, sample_number))
    min_vec = c(0, max_vec[-length(max_vec)])
    segs = mutate(segs,
                  y_min = min_vec[match(ID, factor(unique(segs$ID)))],
                  y_max = max_vec[match(ID, factor(unique(segs$ID)))])
    
    y_labs = distinct(segs, ID, .keep_all = T) %>% 
        mutate(pos = (y_max-y_min)/2 + y_min)
    
    # Plot 
    seg_plot = ggplot(segs, aes(xmin = loc.start, xmax = loc.end, ymin = y_min, ymax = y_max, fill = seg.mean)) +
        geom_rect() +
        scale_fill_gradient2(low = colors[1], mid = colors[2], high = colors[3], guide = 'colourbar',
                             'Log ratio', breaks = legend_breaks, labels = round(legend_breaks)) +
        scale_y_continuous(expand = c(0,0), breaks = y_labs$pos, labels = y_labs$ID) +
        facet_grid(.~chrom, space = 'free_x', scales = 'free_x', switch = 'x') +
        theme(
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.spacing.x = unit(-.5, 'lines'),
            panel.spacing.y = unit(0, 'lines'),
            strip.background = element_rect(fill = 'white'),
            panel.background = element_rect(fill = 'white')
        ) 
    if (return_object == TRUE) {
        seg_plot
    } else {
        print(seg_plot)
    }
}
