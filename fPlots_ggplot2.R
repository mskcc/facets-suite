require(bit64)
require(Cairo)
require(ggplot2)
require(grid)
require(gridExtra)

### Get path to repo
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

copy.number.log.ratio = function(out, fit, load.genome=FALSE, gene.pos=NULL, col.1="#0080FF", col.2="#4CC4FF", sample.num=NULL, lend='butt', theme='bw'){
  
  mat = out$jointseg
  mat = subset(mat, chrom < 23)
  mat = get.cumulative.chr.maploc(mat, load.genome)
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
    geom_point(aes(y=cnlr,x=chr.maploc), colour=c(col.1, col.2)[col.rep], size=.4) +
    scale_x_continuous(breaks=mid, labels=names(mid)) +
    xlab('') +
    ylim(-3,3) +
    ylab('Copy number log ratio') +
    geom_hline(yintercept = dipLogR, color = 'sandybrown', size = .8) +
    geom_segment(data=cncf,aes(x=my.starts$chr.maploc, xend=my.ends$chr.maploc, y=my.starts$cnlr.median, yend=my.ends$cnlr.median), col='red3', size=1, lineend=lend)
  
  panel.grid.col='white'; grid.width = .5
  if(theme=='bw'){panel.grid.col='grey'; grid.width = .2; cnlr = cnlr + theme_bw()}
  
  cnlr = cnlr + theme(axis.text.x  = element_text(angle=90, vjust=0, size=8),
                      axis.text.y = element_text(angle=90, vjust=0, size=8),
                      text = element_text(size=10),
                      panel.grid.minor.x=element_line(colour=panel.grid.col, size=grid.width),
                      panel.grid.major.x=element_line(colour=panel.grid.col, size=0),
                      plot.margin = unit(c(0,1,0,0), 'lines'))  
  cnlr
}

var.allele.log.odds.ratio = function(out, fit, load.genome=FALSE, gene.pos=NULL, col.1="#0080FF", col.2="#4CC4FF", sample.num=NULL, lend='butt', theme='bw'){
  
  mat = out$jointseg
  mat = subset(mat, chrom < 23)
  mat = get.cumulative.chr.maploc(mat, load.genome)
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
    geom_point(aes(y=valor,x=chr.maploc), colour=c(col.1,col.2)[col.rep], size=.4) +
    scale_x_continuous(breaks=mid, labels=names(mid)) +
    xlab('') +
    ylim(-4,4) +
    ylab('Variant allele log odds ratio') +
    geom_segment(data=cncf, aes(x=my.starts$chr.maploc, xend=my.ends$chr.maploc, yend=my.ends$mafR, y=my.starts$mafR), col='red3', size=1, lineend=lend) +
    geom_segment(data=cncf, aes(x=my.starts$chr.maploc, xend=my.ends$chr.maploc, yend=-my.ends$mafR, y=-my.starts$mafR), col='red3', size=1, lineend=lend)
  
  panel.grid.col='white'; grid.width = .5
  if(theme=='bw'){panel.grid.col='grey'; grid.width = .2; valor = valor + theme_bw()}
  
  valor = valor + theme(axis.text.x = element_text(angle=90, vjust=0, size=8),
                        axis.text.y = element_text(angle=90, vjust=0, size=8),
                        text = element_text(size=10),
                        panel.grid.minor.x=element_line(colour=panel.grid.col, size=grid.width),
                        panel.grid.major.x=element_line(colour=panel.grid.col, size=0),
                        plot.margin = unit(c(0,1,0,0), 'lines'))
  valor
}

cellular.fraction = function(out, fit, method=c('cncf', 'em'), load.genome=FALSE, gene.pos=NULL, main='', lend='butt', theme='bw'){
  
  mat = out$jointseg
  mat = subset(mat, chrom < 23)
  mat = get.cumulative.chr.maploc(mat, load.genome)
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
    geom_segment(data=cncf, aes(x=my.starts$chr.maploc, xend=my.ends$chr.maploc, yend=my.ends$cf, y=my.starts$cf), col='black', size=1, lineend=lend) +
    scale_x_continuous(breaks=mid, labels=names(mid)) +
    xlab('') +
    ylim(0,1) +
    ylab(my.ylab)
  
  panel.grid.col='white'; grid.width = .5
  if(theme=='bw'){panel.grid.col='grey'; grid.width = .2; cf = cf + theme_bw()}
  
  cf = cf + theme(axis.text.x  = element_text(angle=90, vjust=0, size=8),
                  axis.text.y = element_text(angle=90, vjust=0, size=8),
                  text = element_text(size=10),
                  panel.grid.minor.x=element_line(colour=panel.grid.col, size=grid.width),
                  panel.grid.major.x=element_line(colour=panel.grid.col, size=0),
                  plot.margin = unit(c(0,1,0,0), 'lines'))
  cf
}

integer.copy.number = function(out, fit, method=c('cncf', 'em'), load.genome=FALSE, gene.pos=NULL, main='', lend='butt', theme='bw'){
  
  mat = out$jointseg
  mat = subset(mat, chrom < 23)
  mat = get.cumulative.chr.maploc(mat, load.genome)
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
    geom_segment(data=cncf, aes(x=my.tcn.starts$chr.maploc, xend=my.tcn.ends$chr.maploc, y=my.tcn.starts$tcn_, yend=my.tcn.ends$tcn_), col='black', size=1,lineend=lend) +
    geom_segment(data=cncf, aes(x=my.lcn.starts$chr.maploc, xend=my.lcn.ends$chr.maploc, y=my.lcn.starts$lcn_, yend=my.lcn.ends$lcn_), col='red', size=1, lineend=lend) +
    scale_y_continuous(breaks=c(0:5, 5 + (1:35)/3), labels=0:40,limits = c(0, NA)) +
    scale_x_continuous(breaks=mid, labels=names(mid)) +
    ylab(my.ylab) +
    xlab('')
  
  panel.grid.col='white'; grid.width = .5
  if(theme=='bw'){panel.grid.col='grey'; grid.width = .2; icn = icn + theme_bw()}
  
  icn = icn + theme(axis.text.x  = element_text(angle=90, vjust=0, size=8),
                    axis.text.y = element_text(angle=90, vjust=0, size=8),
                    text = element_text(size=10),
                    panel.grid.minor.x=element_line(colour=panel.grid.col, size=grid.width),
                    panel.grid.major.x=element_line(colour=panel.grid.col, size=0),
                    plot.margin = unit(c(0,1,0,0), 'lines'))
  
  icn
}

clonal.cluster = function(out, fit, method='em', load.genome=FALSE, gene.pos=NULL, main='', theme='bw') {
  
  mat = out$jointseg
  mat = subset(mat, chrom < 23)
  mat = get.cumulative.chr.maploc(mat, load.genome)
  mid = mat$mid
  mat = mat$mat
  
  cncf = fit$cncf
  cncf = subset(cncf, chrom < 23)
  
  if(method == 'em'){tcnscaled = cncf$tcn.em; tcnscaled[cncf$tcn.em > 5 & !is.na(cncf$tcn.em)] = (5 + (tcnscaled[cncf$tcn.em > 5 & !is.na(cncf$tcn.em)] - 5)/3)}
  if(method == 'cncf'){tcnscaled = cncf$tcn; tcnscaled[cncf$tcn > 5 & !is.na(cncf$tcn)] = (5 + (tcnscaled[cncf$tcn > 5 & !is.na(cncf$tcn)] - 5)/3)}
  tcn_ = rep(tcnscaled, cncf$num.mark)
  
  if(method == 'em'){lcn_ = rep(cncf$lcn.em, cncf$num.mark); my.ylab='CF EM'}
  if(method == 'cncf'){lcn_ = rep(cncf$lcn, cncf$num.mark); my.ylab='Integer copy number (CNCF)'}
  
  mat = cbind(mat, cbind(tcn_, lcn_))
  starts = cumsum(c(1,cncf$num.mark))[1:length(cncf$num.mark)]
  ends = cumsum(c(cncf$num.mark))
  my.tcn.starts = mat[starts,c('chr.maploc','tcn_')]
  my.tcn.ends = mat[ends,c('chr.maploc','tcn_')]
  my.lcn.starts = mat[starts,c('chr.maploc','lcn_')]
  my.lcn.ends = mat[ends,c('chr.maploc','lcn_')]
  
  ccl = ggplot(mat, environment = environment())
  if(!is.null(gene.pos)){
    ccl = ccl + geom_vline(xintercept=gene.pos, color='palevioletred1')
  }
  
  ccf.col = c(colorRampPalette(c("white", "steelblue"))(10),"bisque2")[round(10*cncf$cf.em+0.501)]

  ccl = ccl +
    geom_rect(data=cncf, aes(xmin=my.tcn.starts$chr.maploc, xmax=my.tcn.ends$chr.maploc, ymax=1, ymin=0), fill=ccf.col, col = 'white', size=0) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(breaks=mid, labels=names(mid)) +
    ylab(my.ylab) +
    xlab('')
  
  panel.grid.col='white'; grid.width = .5
  if(theme=='bw'){panel.grid.col='grey'; grid.width = .2; ccl = ccl + theme_bw()}
  
  ccl = ccl + theme(axis.text.x  = element_text(angle=90, vjust=0, size=8),
                    axis.text.y = element_text(angle=90, vjust=0, size=8, color = 'white'),
                    axis.ticks.y = element_line(color = 'white'),
                    text = element_text(size=10),
                    panel.grid.major.y = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    panel.grid.major.x= element_blank(),
                    panel.grid.minor.x= element_blank(),
                    plot.margin = unit(c(0,1,0,0), 'lines'))
  
  ccl
}

get.cumulative.chr.maploc = function(mat, load.genome=FALSE){

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


get.gene.pos = function(hugo.symbol,my.path=paste0(getSDIR(),'/Homo_sapiens.GRCh37.75.canonical_exons.bed'),load.genome=FALSE){

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

  c(mid.point, chrom@values)
}


#Standard facets output plot
plot.facets.all.output = function(out, fit, w=850, h=1100, type='png', load.genome=FALSE, main='', plotname='test',
  gene.name=NULL, lend='butt', em.plot = FALSE) {

  if(!is.null(gene.name)){gene.pos = get.gene.pos(gene.name)}
  if(is.null(gene.name)){gene.pos = NULL}

  if (em.plot == F) {
    cnlr = copy.number.log.ratio(out, fit, gene.pos=gene.pos, lend=lend)
    valor = var.allele.log.odds.ratio(out, fit, gene.pos=gene.pos, lend=lend)
    cfem = cellular.fraction(out, fit, method='em', gene.pos=gene.pos, lend=lend)
    cfcncf = cellular.fraction(out, fit, method='cncf', gene.pos=gene.pos, lend=lend)
    icnem = integer.copy.number(out, fit, method='em', gene.pos=gene.pos, lend=lend)
    icncncf = integer.copy.number(out, fit, method='cncf', gene.pos=gene.pos, lend=lend)  
    all.plots = list(cnlr, valor, cfem, icnem, cfcncf, icncncf)
    plot.no = 6
    plot.h = rep(1, 6)
  } else {
    cnlr = copy.number.log.ratio(out, fit, gene.pos=gene.pos, lend=lend)
    valor = var.allele.log.odds.ratio(out, fit, gene.pos=gene.pos, lend=lend)
    icnem = integer.copy.number(out, fit, method='em', gene.pos=gene.pos, lend=lend)
    cclem = clonal.cluster(out, fit, method = 'em')
    all.plots = list(cnlr, valor, icnem, cclem)
    plot.no = 4
    plot.h = c(1,1,1,.25)
    h = (3/5)*h
  }

  if(type == 'pdf'){plotname = paste(plotname, '.pdf', sep=''); CairoPDF(width = 8.854167, height=11.458333, file=plotname)}
  if(type == 'png'){plotname = paste(plotname, '.png', sep=''); CairoPNG(width = w, height=h, file=plotname, units='px')}

  if(main != ''){grid.arrange(grobs = all.plots,
                              ncol=1,
                              nrow=plot.no,
                              heights = plot.h,
                              top=textGrob(main))}

  if(main == ''){grid.arrange(grobs = all.plots,
                              ncol=1,
                              heights = plot.h,
                              nrow=plot.no)}
  dev.off()

}

#Need to add this functionality so it can be callled by the wrapper, doFacets.R etc.
close.up = function(out, fit, chrom.range=NULL, method=NA, gene.name=NULL, lend='butt', bed.path=NULL){

  if (!is.null(bed.path)) { gene.info = get.gene.pos(gene.name, my.path = bed.path)
  } else { gene.info = get.gene.pos(gene.name) }

  if (!is.null(gene.name)) gene.pos = gene.info[1]
  if (is.null(gene.name)) gene.pos = NULL
  if (is.null(chrom.range)) chrom.range=gene.info[2]
  
  out$out = out$out[out$out$chrom %in% chrom.range,]
  out$jointseg = out$jointseg[out$jointseg$chrom %in% chrom.range,]
  out$IGV = out$IGV[out$IGV$chrom %in% chrom.range,]
  fit$cncf = fit$cncf[fit$cncf$chrom %in% chrom.range,]

  cnlr = copy.number.log.ratio(out, fit, gene.pos=gene.pos, lend=lend)
  valor = var.allele.log.odds.ratio(out, fit, gene.pos=gene.pos, lend=lend)

  output_list <- list(cnlr=cnlr,valor=valor)
  if(method == 'em' | is.na(method)){
    cfem = cellular.fraction(out, fit, method='em', gene.pos=gene.pos, lend=lend)
    icnem = integer.copy.number(out, fit, method='em', gene.pos=gene.pos, lend=lend)
    output_list <- c(output_list, list(cfem=cfem, icnem=icnem))
  }
  if(method == 'cncf' | is.na(method)){
    cfcncf = cellular.fraction(out, fit, method='cncf', gene.pos=gene.pos, lend=lend)
    icncncf = integer.copy.number(out, fit, method='cncf', gene.pos=gene.pos, lend=lend)
    output_list <- c(output_list, list(cfcncf=cfcncf, icncncf=icncncf))
  }
  output_list
}

############################################################################################################################################
############################################################################################################################################

#Example Plot
akt1.close.ups = function(chrom.range = 13:15, gene.name ='AKT1', w=13, h=8, plotname='proj_5513_wxs.pdf',type='pdf',method='cncf'){

  load('~/work//AKT1_UCEC//my_r_003//s_TS01_T/s_TS01_T__s_TS01_N/facets_p300c100/s_TS01_T__s_TS01_N_hisens.Rdata')
  ts01 = close.up(out, fit, chrom.range=chrom.range, method=method, gene.name=gene.name)

  load('~/work//AKT1_UCEC//my_r_003//s_TS02_T/s_TS02_T__s_TS02_N/facets_p300c100/s_TS02_T__s_TS02_N_hisens.Rdata')
  ts02 = close.up(out, fit, chrom.range=chrom.range, method=method, gene.name=gene.name)

  load('~/work//AKT1_UCEC//my_r_003//s_TS03_T/s_TS03_T__s_TS03_N/facets_p300c100/s_TS03_T__s_TS03_N_hisens.Rdata')
  ts03 = close.up(out, fit, chrom.range=chrom.range, method=method, gene.name=gene.name)

  load('~/work//AKT1_UCEC//my_r_003//s_TS04_T/s_TS04_T__s_TS04_N/facets_p300c100/s_TS04_T__s_TS04_N_hisens.Rdata')
  ts04 = close.up(out, fit, chrom.range=chrom.range, method=method, gene.name=gene.name)

  load('~/work//AKT1_UCEC//my_r_003//s_TS05_T/s_TS05_T__s_TS05_N/facets_p300c100/s_TS05_T__s_TS05_N_hisens.Rdata')
  ts05 = close.up(out, fit, chrom.range=chrom.range, method=method, gene.name=gene.name)

  layout = matrix(1:15, nrow = 3)
  if(type == 'pdf'){CairoPDF(width = w, height=h, file=plotname)}
  if(type == 'png'){CairoPNG(width = w, height=h, file=plotname, units='px')}
  grid.arrange(ts01$cnlr, ts02$cnlr, ts03$cnlr, ts04$cnlr, ts05$cnlr,
               ts01$valor, ts02$valor, ts03$valor,  ts04$valor, ts05$valor,
               ts01$icncncf, ts02$icncncf, ts03$icncncf, ts04$icncncf, ts05$icncncf,
               ncol=5, nrow=3)
  dev.off()
}


#Example Plots
akt1.wxs = function(){

  load('~/work//AKT1_UCEC//my_r_003//s_TS01_T/s_TS01_T__s_TS01_N/facets_p300c100/s_TS01_T__s_TS01_N.Rdata')
  #load('/ifs/work/taylorlab/donoghum/AKT1_UCEC//my_r_003//s_TS01_T/s_TS01_T__s_TS01_N/facets_p300c100/s_TS01_T__s_TS01_N_hisens.Rdata')
  plot.facets.all.output(out, fit, type='pdf', main='TS01 | cval: 100', plotname='TS01', gene.name='AKT1')

  load('~/work//AKT1_UCEC//my_r_003//s_TS02_T/s_TS02_T__s_TS02_N/facets_p300c100/s_TS02_T__s_TS02_N.Rdata')
  #load('/ifs/work/taylorlab/donoghum/AKT1_UCEC//my_r_003//s_TS02_T/s_TS02_T__s_TS02_N/facets_p300c100/s_TS02_T__s_TS02_N.Rdata')
  plot.facets.all.output(out, fit, type='pdf', main='TS02 | cval: 100', plotname='TS02', gene.name='AKT1')

  load('~/work//AKT1_UCEC//my_r_003//s_TS03_T/s_TS03_T__s_TS03_N/facets_p300c100/s_TS03_T__s_TS03_N.Rdata')
  #load('/ifs/work/taylorlab/donoghum/AKT1_UCEC//my_r_003//s_TS03_T/s_TS03_T__s_TS03_N/facets_p300c100/s_TS03_T__s_TS03_N.Rdata')
  plot.facets.all.output(out, fit, type='pdf', main='TS03 | cval: 100', plotname='TS03', gene.name='AKT1')

  load('~/work//AKT1_UCEC//my_r_003//s_TS04_T/s_TS04_T__s_TS04_N/facets_p300c100/s_TS04_T__s_TS04_N.Rdata')
  #load('/ifs/work/taylorlab/donoghum/AKT1_UCEC//my_r_003//s_TS04_T/s_TS04_T__s_TS04_N/facets_p300c100/s_TS04_T__s_TS04_N.Rdata')
  plot.facets.all.output(out, fit, type='pdf', main='TS04 | cval: 100', plotname='TS04', gene.name='AKT1')

  load('~/work//AKT1_UCEC//my_r_003//s_TS05_T/s_TS05_T__s_TS05_N/facets_p300c100/s_TS05_T__s_TS05_N.Rdata')
  #load('/ifs/work/taylorlab/donoghum/AKT1_UCEC//my_r_003//s_TS05_T/s_TS05_T__s_TS05_N/facets_p300c100/s_TS05_T__s_TS05_N.Rdata')
  plot.facets.all.output(out, fit, type='pdf', main='TS05 | cval: 100', plotname='TS05', gene.name='AKT1')
}
