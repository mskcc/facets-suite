#!/opt/common/CentOS_6-dev/R/R-3.1.3/bin/Rscript

library(argparse)
library(data.table)
library(ggplot2)
old <- theme_set(theme_bw(30))
library(gtable)
library(gridExtra)
library(grid)

getSDIR <- function(){
  args=commandArgs(trailing=F)
  TAG="--file="
  path_idx=grep(TAG,args)
  SDIR=dirname(substr(args[path_idx],nchar(TAG)+1,nchar(args[path_idx])))
  if(length(SDIR)==0) {
    return(getwd())
  } else {
    return(SDIR)
  }
}

#############################################
### definition of copy number calls in WGD
FACETS_CALL_table <- fread(paste0(getSDIR(), "/FACETS_CALL_table.txt"))
setkey(FACETS_CALL_table, WGD, mcn, lcn)
### lowest value of tcn for AMP
AMP_thresh_tcn <- 6
#############################################



af_grid <- function(maf,
                    output_pdf,
                    log="-",
                    threshold=4,
                    type="ccf",
                    purity=NA,
                    nbins=30,
                    label_WGD = "no WGD"
                    ){
  ### maf must be annotated with tcn and lcn from facets
  maf[,mcn:=as.integer(as.character(tcn))-as.integer(as.character(lcn))]
  mcn_threshold_levels <- c(paste(1:threshold-1), paste0(threshold, "+"))
  maf[,mcn_threshold:=factor(levels = mcn_threshold_levels,
                             ifelse(mcn>=threshold, paste0(threshold, "+"), mcn))]
  maf[,allele_freq:=as.numeric(t_alt_count) / (as.numeric(t_alt_count) + as.numeric(t_ref_count))]
  maf[,allele_freq_bin:=as.numeric(plyr::round_any(allele_freq, 1/nbins))]
  maf[,lcn:=factor(lcn, levels=0:threshold-1)]
  maf <- maf[!is.na(tcn) & !is.na(lcn)]
  maf[, WGD := label_WGD]

  if(is.na(purity))
    purity <- as.numeric(maf[1,list(purity)])
  Tumor_Sample_Barcode <- maf[1,list(Tumor_Sample_Barcode)]

  FACETS_CALL_table <- FACETS_CALL_table[as.integer(mcn) <= 4]
  FACETS_CALL_table[,mcn_threshold:=factor(levels = mcn_threshold_levels,
                                           ifelse(mcn>=threshold, paste0(threshold, "+"), mcn))]
  FACETS_CALL_table <- FACETS_CALL_table[WGD == label_WGD]
  FACETS_CALL_table[, reasonable_icn_status := FALSE]
  FACETS_CALL_table[, sample_purity := purity]
  setkey(FACETS_CALL_table, WGD, mcn_threshold, lcn)

  bin_counts <- maf[, .N, by = list(Tumor_Sample_Barcode, mcn_threshold, lcn, allele_freq_bin)]
  bin_counts[, reasonable_icn_status := ! ( (mcn_threshold == "0" & lcn == "0") | as.integer(as.character(lcn)) > as.integer(as.character(mcn_threshold)) )]
  maxN <<- max(bin_counts$N)
  plot_title <- paste0(Tumor_Sample_Barcode, "   purity: ", round(100*purity), "%")


  p <- ggplot(bin_counts,
              aes(allele_freq_bin,
                  N,
                  fill = reasonable_icn_status)) +
    geom_bar(stat="identity") +
    facet_grid(mcn_threshold ~ lcn, drop=F) +
    scale_x_continuous(limits=c(0,1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    scale_y_continuous(limits=c(0,1.1*maxN)) +
    scale_fill_manual(values=c("red", "black")) +
    ggtitle(plot_title) +
    xlab("Allele Fraction") +
    geom_vline(data=FACETS_CALL_table, aes(xintercept=is.null(FACETS_CALL) + sample_purity), col="blue") +
    geom_vline(data=FACETS_CALL_table, aes(xintercept=is.null(FACETS_CALL) + sample_purity/2), col="blue") +
    #    geom_vline(data=FACETS_CALL_table, aes(xintercept=is.null(FACETS_CALL) + purity/4), col="blue") +
    theme(panel.grid.minor.y = element_blank(),
          legend.position = "none")

  ### add FACETS_CALL description of copy number status
  p <- p + geom_text(data=FACETS_CALL_table, aes(x=0.5, y=1.05 * maxN, label=FACETS_CALL))
  #p

  ### add word "purity" next to the blue line
  purity_annotation_text <- data.table(mcn_threshold = factor(paste0(threshold, "+"), levels = mcn_threshold_levels),
                                       lcn = factor("0", levels=0:threshold-1),
                                       allele_freq_bin = purity,
                                       N = 0,
                                       reasonable_icn_status = FALSE,
                                       lab = "purity")
  p <- p + geom_text(data = purity_annotation_text,
                     aes(label = lab),
                     hjust = ifelse(purity > 0.5, 1.1, -0.1),
                     vjust = -0.1,
                     col = "blue",
                     size = 10)

  ### add Major/Minor Copy Number facet_grid labels
  z <- ggplotGrob(p)
  # add label for right strip
  z <- gtable_add_cols(z,z$widths[[7]])
  z <- gtable_add_grob(z,
                       list(rectGrob(gp = gpar(col = NA, fill = gray(0.5))),
                            textGrob("Major Copy Number", rot = -90, gp = gpar(col = gray(1)))),
                       12, 12, 4, 12, name = paste(runif(2)))

  # add label for top strip
  z <- gtable_add_rows(z, z$heights[[3]], 2)
  z <- gtable_add_grob(z,
                       list(rectGrob(gp = gpar(col = NA, fill = gray(0.5))),
                            textGrob("Minor Copy Number", gp = gpar(col = gray(1)))),
                       3, 4, 3, 10, name = paste(runif(2)))

  # add gaps between facet titles and facet labels
  z <- gtable_add_cols(z, unit(1/8, "line"), 11)
  z <- gtable_add_rows(z, unit(1/8, "line"), 3)

  # draw it
  # grid.newpage()
  z
  # dev.print(dev = pdf, file = output_pdf, width = 14, height = 12)
}

if(!interactive()){

  parser = ArgumentParser()
  parser$add_argument('-m', '--maffile', type='character', help='FACETS-annotated maf to check')
  parser$add_argument('-o', '--outfile', type='character', help='Output filename')
  args=parser$parse_args()

  maf = suppressWarnings(fread(args$maffile))
  outfile = args$outfile
  print(length(unique(maf$Tumor_Sample_Barcode)))
  pdf(outfile, width = 14, height = 12)
   lapply(unlist(unique(maf$Tumor_Sample_Barcode)),
          function(bc, outfile){
            print(grid.arrange(af_grid(maf[Tumor_Sample_Barcode == bc], outfile)))
          },
          outfile = args$outfile)
  # lapply(1:2, function(y){print(grid.arrange(rectGrob(), qplot(y=y)))})
  dev.off()
}
