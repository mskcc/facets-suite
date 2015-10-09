require(Cairo)
require(ggplot2)
require(grid)
require(gridExtra)


copy.number.log.ratio = function(out, fit, load.genome=FALSE, gene.pos=NULL, col.1="#0080FF", col.2="#4CC4FF", sample.num=NULL){
  
  mat = out$jointseg
  mat = subset(mat, chrom < 23)
  mat = get.cumlative.chr.maploc(mat, load.genome)
  mid = mat$mid
  mat = mat$mat
  
  cncf = fit$cncf
  cncf = subset(cncf, chrom < 23)
  dipLogR = out$dipLogR
  
  cnlr.median = rep(cncf$cnlr.median, cncf$num.mark)
  mat = cbind(mat, cnlr.median)
  
  starts = cumsum(c(1,cncf$num.mark))[1:length(cncf$num.mark)]
  ends = cumsum(c(cncf$num.mark))
  my.starts = mat[starts,c('chr.maploc','cnlr.median')]
  my.ends = mat[ends,c('chr.maploc','cnlr.median')]
  
  if(is.null(sample.num)){subset_ = 1:nrow(mat)}
  
  if(is.null(sample.num) == FALSE){
    if(sample.num >= nrow(mat)){subset_ = 1:nrow(mat)}
    if(sample.num < nrow(mat)){subset_ = sort(sample(1:nrow(mat), sample.num, replace=FALSE))}
  }
  
  mat = mat[subset_,]
  col.rep = 1 + rep(mat$chrom - 2 * floor(mat$chrom/2))
  
  cnlr = ggplot(mat, environment = environment()) 
  if(!is.null(gene.pos)){
    cnlr = cnlr + geom_vline(xintercept=gene.pos, color='palevioletred1')
  }
  
  cnlr = cnlr + 
    geom_point(aes(y=cnlr,x=chr.maploc), colour=c(col.1, col.2)[col.rep], size=.8) + 
    scale_x_continuous(breaks=mid, labels=names(mid)) + 
    xlab('') +
    ylim(-3,3) +
    ylab('Log-Ratio') +
    geom_hline(yintercept = dipLogR, color = 'sandybrown', size = .8) + 
    geom_segment(data=cncf,aes(x=my.starts$chr.maploc, xend=my.ends$chr.maploc, y=my.starts$cnlr.median, yend=my.ends$cnlr.median), col='red3', size=1) +
    theme(axis.text.x  = element_text(angle=90, vjust=0, size=8),
          axis.text.y = element_text(angle=90, vjust=0, size=8),
          text = element_text(size=10))
  cnlr
}

var.allele.log.odds.ratio = function(out, fit, load.genome=FALSE, gene.pos=NULL, col.1="#0080FF", col.2="#4CC4FF", sample.num=NULL){
  
  mat = out$jointseg
  mat = subset(mat, chrom < 23)
  mat = get.cumlative.chr.maploc(mat, load.genome)
  mid = mat$mid
  mat = mat$mat
  
  cncf = fit$cncf
  cncf = subset(cncf, chrom < 23)
  
  mafR = rep(sqrt(abs(cncf$mafR)), cncf$num.mark)
  mat = cbind(mat, mafR)
  
  starts = cumsum(c(1,cncf$num.mark))[1:length(cncf$num.mark)]
  ends = cumsum(c(cncf$num.mark))
  my.starts = mat[starts,c('chr.maploc','mafR')]
  my.ends = mat[ends,c('chr.maploc','mafR')]
  
  if(is.null(sample.num)){subset_ = 1:nrow(mat)}
  
  if(is.null(sample.num) == FALSE){
    if(sample.num >= nrow(mat)){subset_ = 1:nrow(mat)}
    if(sample.num < nrow(mat)){subset_ = sort(sample(1:nrow(mat), sample.num, replace=FALSE))}
  }  
  
  mat = mat[subset_,]
  col.rep = 1 + rep(mat$chrom - 2 * floor(mat$chrom/2))
  
  valor = ggplot(mat, environment = environment()) 
  if(!is.null(gene.pos)){
    valor = valor + geom_vline(xintercept=gene.pos, color='palevioletred1')
  }
  
  valor = valor +
    geom_point(aes(y=valor,x=chr.maploc), colour=c(col.1,col.2)[col.rep], size=.8) +
    scale_x_continuous(breaks=mid, labels=names(mid)) +  
    xlab('') +
    ylim(-4,4) +
    ylab('Log-Odds-Ratio') +
    geom_segment(data=cncf, aes(x=my.starts$chr.maploc, xend=my.ends$chr.maploc, yend=my.ends$mafR, y=my.starts$mafR), col='red3', size=1) +
    geom_segment(data=cncf, aes(x=my.starts$chr.maploc, xend=my.ends$chr.maploc, yend=-my.ends$mafR, y=-my.starts$mafR), col='red3', size=1) +
    theme(axis.text.x = element_text(angle=90, vjust=0, size=8),
          axis.text.y = element_text(angle=90, vjust=0, size=8),
          text = element_text(size=10))
  valor
}

cellular.fraction = function(out, fit, method=c('cncf', 'em'), load.genome=FALSE, gene.pos=NULL, main=''){
  
  mat = out$jointseg
  mat = subset(mat, chrom < 23)
  mat = get.cumlative.chr.maploc(mat, load.genome)
  mid = mat$mid
  mat = mat$mat
  
  cncf = fit$cncf
  cncf = subset(cncf, chrom < 23)
  
  if(method == 'em'){cncf$cf.em[is.na(cncf$cf.em)] = -1; cf = rep(cncf$cf.em, cncf$num.mark); my.ylab='Cellular fraction (EM)'}  
  if(method == 'cncf'){cncf$cf[is.na(cncf$cf)] = -1; cf = rep(cncf$cf, cncf$num.mark); my.ylab='Cellular fraction (CNCF)'}
  
  mat = cbind(mat, cf)
  starts = cumsum(c(1,cncf$num.mark))[1:length(cncf$num.mark)]
  ends = cumsum(c(cncf$num.mark))
  my.starts = mat[starts,c('chr.maploc','cf')]
  my.ends = mat[ends,c('chr.maploc','cf')]
  
  cf = ggplot(mat, environment = environment()) 
  if(!is.null(gene.pos)){
    cf = cf + geom_vline(xintercept=gene.pos, color='palevioletred1')
  }
  
  cf = cf +
    geom_segment(data=cncf, aes(x=my.starts$chr.maploc, xend=my.ends$chr.maploc, yend=my.ends$cf, y=my.starts$cf), col='black', size=1) +
    scale_x_continuous(breaks=mid, labels=names(mid)) +
    xlab('') +
    ylim(0,1) +
    ylab(my.ylab) +
    theme(axis.text.x  = element_text(angle=90, vjust=0, size=8),
          axis.text.y = element_text(angle=90, vjust=0, size=8),
          text = element_text(size=10))
  cf
}

integer.copy.number = function(out, fit, method=c('cncf', 'em'), load.genome=FALSE, gene.pos=NULL, main=''){
  
  mat = out$jointseg
  mat = subset(mat, chrom < 23)
  mat = get.cumlative.chr.maploc(mat, load.genome)
  mid = mat$mid
  mat = mat$mat
  
  cncf = fit$cncf
  cncf = subset(cncf, chrom < 23)
  
  if(method == 'em'){tcnscaled = cncf$tcn.em; tcnscaled[cncf$tcn.em > 5 & !is.na(cncf$tcn.em)] = (5 + (tcnscaled[cncf$tcn.em > 5 & !is.na(cncf$tcn.em)] - 5)/3)}
  if(method == 'cncf'){tcnscaled = cncf$tcn; tcnscaled[cncf$tcn > 5 & !is.na(cncf$tcn)] = (5 + (tcnscaled[cncf$tcn > 5 & !is.na(cncf$tcn)] - 5)/3)}
  tcn_ = rep(tcnscaled, cncf$num.mark)
  
  if(method == 'em'){lcn_ = rep(cncf$lcn.em, cncf$num.mark); my.ylab='Integer copy number (EM)'}
  if(method == 'cncf'){lcn_ = rep(cncf$lcn, cncf$num.mark); my.ylab='Integer copy number (CNCF)'}
  
  mat = cbind(mat, cbind(tcn_, lcn_))
  starts = cumsum(c(1,cncf$num.mark))[1:length(cncf$num.mark)]
  ends = cumsum(c(cncf$num.mark))
  my.tcn.starts = mat[starts,c('chr.maploc','tcn_')]
  my.tcn.ends = mat[ends,c('chr.maploc','tcn_')]
  my.lcn.starts = mat[starts,c('chr.maploc','lcn_')]
  my.lcn.ends = mat[ends,c('chr.maploc','lcn_')]
  
  icn = ggplot(mat, environment = environment()) 
  if(!is.null(gene.pos)){
    icn = icn + geom_vline(xintercept=gene.pos, color='palevioletred1')
  }
  
  icn = icn + 
    geom_segment(data=cncf, aes(x=my.tcn.starts$chr.maploc, xend=my.tcn.ends$chr.maploc, y=my.tcn.starts$tcn_, yend=my.tcn.ends$tcn_), col='black', size=1) +
    geom_segment(data=cncf, aes(x=my.lcn.starts$chr.maploc, xend=my.lcn.ends$chr.maploc, y=my.lcn.starts$lcn_, yend=my.lcn.ends$lcn_), col='red', size=1) +
    scale_y_continuous(breaks=c(0:5, 5 + (1:35)/3), labels=0:40,limits = c(0, NA)) + 
    scale_x_continuous(breaks=mid, labels=names(mid)) + 
    ylab(my.ylab) +
    xlab('') +
    theme(axis.text.x  = element_text(angle=90, vjust=0, size=8),
          axis.text.y = element_text(angle=90, vjust=0, size=8),
          text = element_text(size=10))
  icn
}

get.cumlative.chr.maploc = function(mat, load.genome=FALSE){

  if(load.genome){
    require(BSgenome.Hsapiens.UCSC.hg19)
    genome = BSgenome.Hsapiens.UCSC.hg19
    chrom.lengths = seqlengths(genome)[1:22]
  }
  else{
    chrom.lengths = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,135006516, 
                      133851895, 115169878, 107349540, 102531392, 90354753,  81195210,  78077248,  59128983,  63025520, 48129895,  51304566) #hg19
  }
  
  cum.chrom.lengths = cumsum(as.numeric(chrom.lengths))
  mid = cum.chrom.lengths - (chrom.lengths/2)
  names(mid) = 1:22
  
  chr.maploc.to.gen.maploc = function(x){mat[mat$chrom==x,]$maploc + cum.chrom.lengths[x-1]}
  chr.maploc = sapply(2:22,chr.maploc.to.gen.maploc)
  chr.maploc = unlist(chr.maploc)
  chr.maploc = c(mat[mat$chrom==1,]$maploc,chr.maploc)
  mat = cbind(mat,chr.maploc)
  
  list(mat=mat, mid=mid)
}

get.gene.pos = function(hugo.symbol,my.path='~/home/reference_sequences/Homo_sapiens.GRCh37.75.canonical_exons.bed',load.genome=FALSE){
  
  if(load.genome){
    require(BSgenome.Hsapiens.UCSC.hg19)
    genome = BSgenome.Hsapiens.UCSC.hg19
    chrom.lengths = seqlengths(genome)[1:22]
  }
  else{
    chrom.lengths = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,135006516, 
                      133851895, 115169878, 107349540, 102531392, 90354753,  81195210,  78077248,  59128983,  63025520, 48129895,  51304566) #hg19
  }
  
  cum.chrom.lengths = cumsum(as.numeric(chrom.lengths))
  
  require(rtracklayer)
  genes = import.bed(my.path)
  mcols(genes)$name = matrix(unlist(strsplit(mcols(genes)$name,':')),nc=3,byrow=T)[,1]
  gene.start = min(start(genes[which(mcols(genes)$name == hugo.symbol)]))
  gene.end = max(end(genes[which(mcols(genes)$name == hugo.symbol)]))
  mid.point = gene.start + ((gene.end - gene.start)/2) 
  chrom = seqnames(genes[which(mcols(genes)$name == hugo.symbol)])[1]  
  mid.point = cum.chrom.lengths[as.integer(chrom)-1] + mid.point
  
  mid.point
}

#Standard facets output plot
plot.facets.all.output = function(out, fit, w=850, h=1100, type='png', load.genome=FALSE, main='', plotname='test', gene.name=NULL){
  
  if(!is.null(gene.name)){gene.pos = get.gene.pos(gene.name)}
  
  cnlr = copy.number.log.ratio(out, fit, gene.pos=gene.pos)
  valor = var.allele.log.odds.ratio(out, fit, gene.pos=gene.pos)
  cfem = cellular.fraction(out, fit, method='em', gene.pos=gene.pos)
  cfcncf = cellular.fraction(out, fit, method='cncf', gene.pos=gene.pos)
  icnem = integer.copy.number(out, fit, method='em', gene.pos=gene.pos)
  icncncf = integer.copy.number(out, fit, method='cncf', gene.pos=gene.pos)
  
  if(type == 'pdf'){plotname = paste(plotname, '.pdf', sep=''); CairoPDF(width = 8.854167, height=11.458333, file=plotname)}
  if(type == 'png'){plotname = paste(plotname, '.png', sep=''); CairoPNG(width = w, height=h, file=plotname, units='px')}
  
  if(main != ''){grid.arrange(cnlr, valor, cfem, icnem, cfcncf, icncncf,
                                  ncol=1,
                                  nrow=6,
                                  top=textGrob(main))}
  
  if(main == ''){grid.arrange(cnlr, valor, cfem, icnem, cfcncf, icncncf,
                                  ncol=1,
                                  nrow=6)}
  dev.off()
  
}

#Need to add this functionality so it can be callled by the wrapper, doFacets.R etc.
close.up = function(out, fit, chrom.range, method=NULL, gene.name=NULL){
  
  out$out = out$out[out$out$chrom %in% chrom.range,]
  out$jointseg = out$jointseg[out$jointseg$chrom %in% chrom.range,] 
  out$IGV = out$IGV[out$IGV$chrom %in% chrom.range,]
  fit$cncf = fit$cncf[fit$cncf$chrom %in% chrom.range,]
  
  if(!is.null(gene.name)){gene.pos = get.gene.pos(gene.name)}
  
  cnlr = copy.number.log.ratio(out, fit, gene.pos=gene.pos)
  valor = var.allele.log.odds.ratio(out, fit, gene.pos=gene.pos)
  cfem = cellular.fraction(out, fit, method='em', gene.pos=gene.pos)
  cfcncf = cellular.fraction(out, fit, method='cncf', gene.pos=gene.pos)
  icnem = integer.copy.number(out, fit, method='em', gene.pos=gene.pos)
  icncncf = integer.copy.number(out, fit, method='cncf', gene.pos=gene.pos)
  
  list(cnlr=cnlr,valor=valor,cfem=cfem,cfcncf=cfcncf,icnem=icnem,icncncf=icncncf)
}
