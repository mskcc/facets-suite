#!/usr/bin/env Rscript

##########################################################################################
##########################################################################################
# MSKCC CMO
# FACETS QC
# Compare FACETS purity estimate to GMM of VAF distribution in maf file
##########################################################################################
##########################################################################################

'%!in%' <- function(x,y)!('%in%'(x,y))

catverbose <- function(...) {
  cat(format(Sys.time(), "%Y%m%d %H:%M:%S |"), ..., "\n")
}

monotonic <- function(x) {
  all(x == cummax(x))
}

plot_vaf_by_cn_state <- function(maf, sample, purity, wgd=F) {
  maf.tmp <- maf[lcn == mcn]
  gg <- ggplot(maf.tmp, aes(x = VAF)) + 
    geom_histogram(col = "black", fill="#41B6C4", lwd = 1.5, binwidth = 0.02) +
    geom_vline(xintercept = (purity/2), 
               linetype = 2, 
               color = "#FB6A4A") +
    xlim(c(0,1)) +
    facet_grid(lcn ~ mcn) + 
    xlab("Variant Allele Fraction") +
    ylab("Frequency") +
    ggtitle(sample) +
    theme_bw() + 
    theme(plot.title=element_text(size=25, face = "bold"),
          axis.title=element_text(size=20, face = "bold"),
          strip.text.x=element_text(size=20, face = "bold"),
          strip.text.y=element_text(size=20, face = "bold"),
          axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15))
  if (wgd) {
    gg <- gg + geom_vline(xintercept = (purity/4), linetype = 2, color = "#FD8D3C")
  }
  plot(gg)
}

main <- function(maf, facets, plot=F) {
  
  summary <- c()
  samples <- unique(maf$Tumor_Sample_Barcode)
  
  for (s in samples) {
  
    cat('\n')
    catverbose(s)
    s.maf <- maf[Tumor_Sample_Barcode == s & mcn == lcn]
    s.facets <- facets[Tumor_Sample_Barcode == s]$Rdata_filename
    wgd <- F
    flag <- F
    
    catverbose("Loading FACETS Rdata...")
    load(s.facets)
    facets.fit <- as.data.table(fit$cncf)
    dipLogR <- out$dipLogR
    s.purity <- fit$purity
    
    n.bases <- sum(as.numeric(fit$seglen))
    if(n.bases > 0) {
      f_hi_mcn <- sum(as.numeric(fit$seglen[which((facets.fit$tcn - facets.fit$lcn) >= 2)])) / n.bases
    } else {
      f_hi_mcn <- NA
    }
    
    if(!is.na(f_hi_mcn) & f_hi_mcn > 0.5) { # Major copy number >= 2 across 50% of the genome
      wgd <- T
    }
    
    if (plot) {
      pdf(file=paste0(s, "_balanced_vafs.pdf"), width=20, height=14)
      plot_vaf_by_cn_state(s.maf, s, s.purity, wgd)
      dev.off()
    }
    
    # mclust.options(classPlotColors = brewer.pal(n=9, "YlGnBu")[c(2,5,8)])
    catverbose("Fitting Gaussian Mixture Model...")
    
    # Specify extra component for sublconal mutations
    vaf_min <- min(s.maf$VAF)
    if (!wgd) {
      prior.means <- c(vaf_min, s.purity/2)
      s.comp <- 2
    } else {
      prior.means <- c(vaf_min, s.purity/4, s.purity/2)
      s.comp <- 3
    }
    
    dmclust <- densityMclust(s.maf$VAF, 
                             prior=priorControl(mean=prior.means),
                             G=s.comp) 
    
    s.maf[, clust := dmclust$classification]
    if (monotonic(s.maf[order(VAF)]$clust)) {
      
      summary(dmclust, parameters = T)
      
      if (plot) {
        pdf(file=paste0(s, "_facets_vaf_test.pdf"), width=20, height=14)
        par(mfrow=c(3,1))
        plot(dmclust, what = "density", data = s.maf$VAF, xlab = "VAF", breaks=1:100/100)
        plot(dmclust, what = "diagnostic", type="qq")
        plot(dmclust, what = "diagnostic", type="cdf")
        dev.off()
      }
      
      if (dmclust$modelName == "E") {
        dmclust$parameters$variance$sigmasq <- rep(dmclust$parameters$variance$sigmasq[1],
                                                   times=s.comp)
      }
      
      # Use t-tests to allow for clusters w/ small n
      n1 <- sum(dmclust$classification == 1)
      mu1 <- dmclust$parameters$mean[1]
      sd1 <- sqrt(dmclust$parameters$variance$sigmasq[1])
      t1 <- (mu1 - s.purity/2) / (sd1 / sqrt(n1))
      pval1 <- 2*pt(-abs(t1), df = n1-1)
      
      n2 <- sum(dmclust$classification == 2)
      mu2 <- dmclust$parameters$mean[2]
      sd2 <- sqrt(dmclust$parameters$variance$sigmasq[2])
      t2 <- (mu2 - s.purity/2) / (sd2 / sqrt(n2))
      pval2 <- 2*pt(-abs(t2), df = n2-1)
      
      pvals <- c(pval1, pval2)
      if (s.comp == 3) {
        
        n3 <- sum(dmclust$classification == 2)
        mu3 <- dmclust$parameters$mean[2]
        sd3 <- sqrt(dmclust$parameters$variance$sigmasq[2])
        t3 <- (mu2 - s.purity/2) / (sd2 / sqrt(n2))
        pval3 <- 2*pt(-abs(t3), df = n3-1)
        pvals <- c(pvals, pval3)
        
      }
      
        catverbose(paste0("Cluster p-value: ", pvals))
        if(all(pvals < .01/s.comp)) { # Conservative adjustment
          flag <- T
        }
    } else {
      catverbose("Bad model fit")
    }
    
    if (flag) {
      alt.purity <- unname(2*dmclust$parameters$mean[which.max(dmclust$parameters$pro)])
    } else {
      alt.purity <- NA
    }
    
    s.summary <- c(s, s.purity, dipLogR, wgd, alt.purity)
    summary <- rbind(summary, s.summary)
    
  }
  
  rownames(summary) <- NULL
  colnames(summary) <- c('Tumor_Sample_Barcode', 'Purity', 'dipLogR', 'WGD', 'Alt_Purity')
  write.table(summary, file="FACETS_VAF_Summary.txt", quote=F, row.names=F, col.names=T, sep="\t")
  
}


if ( ! interactive() ) {
  
  pkgs = c('argparse', 'data.table', 'mclust', 'RColorBrewer', 'ggplot2')
  tmp <- lapply(pkgs, function (x){
    suppressPackageStartupMessages(require(x, character.only = TRUE))
  })
  rm(tmp)
  options(datatable.showProgress = FALSE)
  
  parser=ArgumentParser()
  parser$add_argument('-m', '--maf', required=T, type='character', help='MAF file with FACETS annotations')
  parser$add_argument('-f', '--facets', required=T, type='character', 
                      help='Tab-delimited mapping of "Tumor_Sample_Barcode" from maf file to FACETS Rdata files')
  parser$add_argument('-p', '--plot', action="store_true", default=F, help='Plot model output')
  args=parser$parse_args()
  
  maf <- fread(args$maf)
  facets <- fread(args$facets)
  plot <- args$plot
  
  if ('mcn' %!in% names(maf)) {
    maf[, mcn := tcn-lcn]
  }
  
  sink('FACETS_FIT_QC.log')
  main(maf, facets, plot)
  sink()
  
}
