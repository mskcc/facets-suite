#!/usr/bin/env Rscript

##########################################################################################
##########################################################################################
# MSKCC CMO
# FACETS QC
##########################################################################################
##########################################################################################

'%!in%' <- function(x,y)!('%in%'(x,y))

catverbose <- function(...){
  cat(format(Sys.time(), "%Y%m%d %H:%M:%S |"), ..., "\n")
}

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

plot_vaf_by_cn_state <- function(maf, wgd=F) {
  
  if ('mcn' %!in% names(maf)){
    maf[, mcn := tcn-lcn]
  }
  
  maf.tmp <- maf[!is.na(mcn) & mcn <= 6]
  if (!('t_var_freq' %in% names(maf.tmp))) maf.tmp[,t_var_freq := t_alt_count/t_depth]
  phi <- unique(maf[!is.na(purity)]$purity)
  if(length(phi) == 0){
    catverbose("No FACETS annotations!")
  } else {
    gg <- ggplot(maf.tmp, aes(x=t_var_freq)) + 
      geom_histogram(col="black", fill="#41B6C4", lwd=1.5, binwidth = 0.02) +
      geom_vline(xintercept=(phi/2), linetype=2, color = "#FB6A4A") +
      facet_grid(lcn ~ mcn) + 
      xlab("Variant Allele Fraction") +
      ylab("Frequency") +
      theme_bw() + 
      theme(plot.title=element_text(size=25, face = "bold"),
            axis.title=element_text(size=20, face = "bold"),
            strip.text.x=element_text(size=20, face = "bold"),
            strip.text.y=element_text(size=20, face = "bold"),
            axis.text.x=element_text(size=15, angle=45, hjust=1),
            axis.text.y=element_text(size=15),
            legend.text=element_text(size=15),
            legend.title=element_text(size=15))
    if(wgd){
      gg <- gg + geom_vline(xintercept=(phi/4), linetype=2, color = "#FD8D3C")
    }
    gg
  }
}

center_igv_file <- function(outfile){
  
  igv.adj <- out$IGV
  
  if(out$dipLogR <= 0){igv.adj$seg.mean = igv.adj$seg.mean + abs(out$dipLogR)}
  if(out$dipLogR > 0){igv.adj$seg.mean = igv.adj$seg.mean - out$dipLogR}
  
  write.table(igv.adj, file = outfile, quote = F, row.names = F, col.names = T,
              sep = "\t")
  
}

facets_qc <- function(maf, facets, igv=F){
  
  summary <- c()
  flagged <- c()
  samples <- unique(maf$Tumor_Sample_Barcode)
  
  for (s in samples){
    
    cat('\n')
    catverbose(s)
    s.maf <- maf[Tumor_Sample_Barcode == s]
    s.facets <- facets[Tumor_Sample_Barcode == s]$Rdata_filename
    alt.fit <- F
  
    load(s.facets)
    catverbose("Loading FACETS Rdata...")
    facets.fit <- as.data.table(fit$cncf)
    
    purity <- fit$purity
    ploidy <- fit$ploidy
    catverbose(paste0("Purity: ", purity, " | ", "Ploidy: ", ploidy))
    
    if(!is.null(out$flags)){
      alt.fit <- T
      catverbose(paste0("Flags: ", out$flags))
    }
    
    if(purity < 0.3 | is.na(purity)){
      catverbose(paste0("Purity < 0.3, use EM"))
    }
    
    ### TODO: report egregious mismatches between EM, CNCF
    
    dipLogR <- out$dipLogR
    if(abs(dipLogR) > 1){
      alt.fit <- T
      catverbose(paste0("dipLogR magnitude > 1"))
    }
    
    # WGD
    wgd <- F
    f_hi_mcn <- sum(as.numeric(fit$seglen[which((facets.fit$tcn - facets.fit$lcn) >= 2)])) / sum(as.numeric(fit$seglen))
    if(f_hi_mcn > 0.5){ # Major copy number >= 2 across 50% of the genome
      wgd <- T
      alt.fit <- T
      catverbose(paste0("Likely WGD"))
    }
    
    # Balanced diploid region in WGD case 
    balance.thresh <- out$mafR.thresh
    if(wgd){
      dipLogR.bal.segs <- facets.fit[cnlr.median.clust == dipLogR & mafR.clust < balance.thresh]
      dipLogR.bal.chrs <- unique(dipLogR.bal.segs$chrom)
      dipLogR.bal.out <- paste(dipLogR.bal.chrs, collapse = "|")
      if(length(dipLogR.bal.chrs > 0)){
        catverbose("Balanced segments at dipLogR in chromosomes:")
        catverbose(dipLogR.bal.chrs)
      }
    }
    
    # LOH
    loh <- sum(as.numeric(fit$seglen[which(facets.fit$lcn == 0)])) / sum(as.numeric(fit$seglen))
    if(loh > 1/2){
      catverbose(paste0("Widespread LOH"))
    }
    
    # Extreme amplifications and homozygous deletions
    n.amps <- nrow(facets.fit[tcn > 10])
    n.homdels <- nrow(facets.fit[tcn == 0])
    
    # Hypersegmentation / complex rearrangements
    chr.break <- F
    n.segments <- table(unique(facets.fit[, .(chrom, tcn, lcn)])[["chrom"]])
    if(any(n.segments >= 5)){ # 5 or more unique copy number states in a given chromosome 
      chr.break <- T
      chr.tmp <- names(which(table(unique(facets.fit[, .(chrom, tcn, lcn)])[["chrom"]]) > 5))
      catverbose("5+ unique CN states on chromosome(s):") 
      catverbose(chr.tmp)
    }
    
    bal.logR <- out$alBalLogR[, 'dipLogR']
    if(dipLogR %!in% bal.logR){
      catverbose("Unbalanced dipLogR")
    }
    
    if(alt.fit){
      alt.diplogR <- as.numeric(bal.logR[which(bal.logR != dipLogR)])
      catverbose(paste0("Alternative dipLogR: ", alt.diplogR))
    }
    
    if (!(all(is.na(s.maf$lcn)))) {# Faceting of plot won't work if this is the case
      ggsave(plot_vaf_by_cn_state(s.maf, wgd), filename = paste0(s, "_vaf_vs_cn.pdf"),
        width = 20, height = 14)
    }
    
    homloss.idx <- which(facets.fit$tcn == 0)
    facets.fit.homloss <- as.data.table(cbind(facets.fit$chrom[homloss.idx], 
                                              fit$start[homloss.idx],
                                              fit$end[homloss.idx])
                                        )
    
    setnames(facets.fit.homloss, c("Chromosome", "Start_Position", "End_Position"))
    setkey(facets.fit.homloss, Chromosome, Start_Position, End_Position)
    facets.fit.homloss <- facets.fit.homloss[Chromosome %in% seq(1,22)]
    
    s.maf.tmp <- s.maf[Chromosome %in% seq(1,22), .(Chromosome, Start_Position, End_Position)]
    s.maf.tmp[, Chromosome := as.numeric(Chromosome)]
    
    homloss.overlaps <- foverlaps(s.maf.tmp,
                                  facets.fit.homloss,
                                  type = "any")[!is.na(Start_Position)]
    if(nrow(homloss.overlaps) > 0){
      catverbose(paste0("Mutations called in homozygous loss segments"))
    }
    
    if(igv){
      center_igv_file(outfile = paste0(s, "_igv.adj.seg"))
    }
    
    s.summary <- c(s, purity, ploidy, dipLogR, f_hi_mcn, wgd, loh, n.amps, n.homdels)
    summary <- rbind(summary, s.summary)
    
    if(!is.null(out$flags)) {
      flags.out <- paste(out$flags, collapse = " | ")
      s.summary <- c(s.summary, flags.out)
      flagged <- rbind(flagged, s.summary) 
    }
  
  }
  
  colnames(summary) <- c('Tumor_Sample_Barcode', 'Purity', 'Ploidy', 'dipLogR', 'f_hi_MCN', 'WGD', 'LOH', 'Amps.n', 'HomDels.n')
  write.table(summary, file="FACETS_QC_summary.txt", quote=F, row.names=F, col.names=T,
              sep="\t")
  
  colnames(flagged) <- c('Tumor_Sample_Barcode', 'Purity', 'Ploidy', 'dipLogR', 'f_hi_MCN', 'WGD', 'LOH', 'Amps.n', 'HomDels.n', 'Flags')
  write.table(flagged, file="FACETS_QC_flagged.txt", quote=F, row.names=F, col.names=T,
              sep="\t")
  

}

if ( ! interactive() ) {
  
  pkgs = c('data.table', 'argparse', 'RColorBrewer', 'ggplot2')
  tmp <- lapply(pkgs, require, quietly=T, character.only = T)
  rm(tmp)
  options(datatable.showProgress = FALSE)
  
  parser=ArgumentParser()
  parser$add_argument('-m', '--maf', required=T, type='character', help='MAF file with FACETS annotations')
  parser$add_argument('-f', '--facets', required=T, type='character', 
                      help='Mapping of "Tumor_Sample_Barcode" from maf and "Rdata_filename" from FACETS (tab-delimited with header)')
  parser$add_argument('-i', '--igv', action='store_true', default=F,
                      help='Output adjusted seg file for IGV')
  args=parser$parse_args()

  maf <- fread(args$maf)
  facets <- fread(args$facets)
  setnames(facets, c("Tumor_Sample_Barcode", "Rdata_filename"))
  igv <- args$igv
  
  sink('FACETS_QC.log')
  facets_qc(maf, facets, igv)
  sink()
  
}


