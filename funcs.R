runEM2=function(d,matched.normal=TRUE) {
	pmat=d$pmat
	jointseg=d$jointseg
	out=d$out
	jointseg=subset(jointseg,!is.na(cnlr)); #remove positions with missing logR
	out=data.frame(out,cnlr.mean=out$cnlr.median); #temp fix. variable name change in naivegcf
	fit=emgcf(jointseg,out,matched.normal=matched.normal,trace=F,maxiter=10,smooth=T,dipLogR=-0.74)
	out3=naivegcf(out,fit$dipLogR)
	combined=data.frame(
		fit$out2,cf.naivegcf=out3$cf,tcn.naivgcf=out3$tcn,lcn.naivegcf=out3$lcn,
		purity=fit$purity,ploidy=fit$ploidy,wgd=ifelse(fit$dipt==4,1,0)
	)
	return(list(pmat=pmat,jointseg=jointseg,out=out,out2=fit$out2,dipLogR=fit$dipLogR,combined=combined))	
}

# To resegment previously segmented samples (by Nick), provided the Rdata file
# exists of the original output. Testing the increase of cval (default=35 for
# exome-wide data) for WGS data.
reSeg=function(orig,cval=35) {
	dmat=counts2logROR(orig$pmat)
	jointseg=segsnps(dmat,cval)
	out=jointsegsummary(jointseg)
	out1=clustersegs(jointseg,out,25,cval)
	jointseg$segclust=rep(NA_real_,nrow(jointseg))
	jointseg$segclust[is.finite(jointseg$cnlr)]=out1$snpclust
	out$segclust=rep(NA_real_,nrow(out))
	out$segclust=out1$segclust
	out$cnlr.median.clust=out1$segclustsummary[out$segclust,"cnlr.median"]
	out$mafR.clust=out1$segclustsummary[out$segclust,"mafR"]
	list(pmat=orig$pmat,jointseg=jointseg,out=out)
}

makeLegacySegObject=function(out,sampID) {
	suppressMessages(require(DNAcopy))
	d=list()
	d$data=out$jointseg[,c("chrom","maploc","cnlr")]
	colnames(d$data)[3]=sampID
	d$output=formatSegmentOutput(out,sampID)
	d$call="facets.seg"
	class(d$data)=c("CNA","data.frame")
	class(d)="DNAcopy"
	d
}

formatSegmentOutput=function(out,sampID) {
	seg=list()
	seg$ID=rep(sampID,nrow(out$out))
	seg$chrom=out$out$chr
	seg$loc.start=rep(NA,length(seg$ID))
	seg$loc.end=seg$loc.start
	seg$num.mark=out$out$num.mark
	seg$seg.mean=out$out$cnlr.median
	for(i in 1:nrow(out$out)) {
#		lims=range(out$jointseg$maploc[(out$jointseg$chrom==out$out$chr[i] & out$jointseg$segs==out$out$seg[i])],na.rm=T)
		lims=range(out$jointseg$maploc[(out$jointseg$chrom==out$out$chrom[i] & out$jointseg$seg==out$out$seg[i])],na.rm=T)
		seg$loc.start[i]=lims[1]
		seg$loc.end[i]=lims[2]
	}
	as.data.frame(seg)
}

plotSample3=function(x,pch=".",cex=2,label="NA") {
	def.par=par(no.readonly = TRUE)
	jseg=x$jointseg
	out=x$out
	chrcol=1+rep(out$chr - 2 * floor(out$chr/2),out$num.mark)
	nn=cumsum(table(x$jointseg$chrom))
	segbdry=cumsum(c(0,out$num.mark))
	segstart=segbdry[-length(segbdry)]
	segend=segbdry[-1]
	layout(matrix(c(1,1,2,2),ncol=1))
	par(mar=c(0.25,3,0.25,1),mgp=c(2,0.7,0),oma=c(3,0,1.25,0))
	colors=c("grey","lightblue","azure4","slateblue")
	plot(jseg$cnlr[!is.na(jseg$cnlr)],pch=pch,cex=cex,col=colors[chrcol],ylab="log-ratio",xaxt="n")
	abline(h=median(jseg$cnlr,na.rm=TRUE),col="green2")
	segments(segstart,out$cnlr.median,segend,out$cnlr.median,lwd=1.75,col=2)
	mtext(side=3,line=0,at=length(jseg$cnlr[!is.na(jseg$cnlr)])/2,paste("Sample",label),cex=0.8)
	#if(length(nn)==23) {
	#	mtext(c(1:22,"X"),side=3,line=0,at=(nn+c(0,nn[-23]))/2,cex=0.65)
	#} else {
	#	mtext(1:22,side=3,line=0,at=(nn+c(0,nn[-22]))/2,cex=0.65)
	#}
	plot(jseg$valor[!is.na(jseg$cnlr)],pch=pch,cex=cex+0.5,col=colors[chrcol],ylab="log-odds-ratio",ylim=c(-5,5),xaxt="n")
	segments(segstart,sqrt(abs(out$mafR)),segend,sqrt(abs(out$mafR)),lwd=1.75,col=2)
	segments(segstart,-sqrt(abs(out$mafR)),segend,-sqrt(abs(out$mafR)),lwd=1.75,col=2)
	#axis(1)
	#mtext(side=1,line=1.75,label,cex=0.8)
	if(length(nn)==23) { 
		mtext(c(1:22,"X"),side=1,line=0,at=(nn+c(0,nn[-23]))/2,cex=0.65)
	} else { 
		mtext(1:22,side=1,line=0,at=(nn+c(0,nn[-22]))/2,cex=0.65)
	}
	par(def.par)
}

plotSampleCNCF2=function(out,fit) {
	orig=out
  mat=out$jointseg
  mat=subset(mat,chrom<23)
  out=subset(out$out,chr<23)
  cncf=fit$cncf
  cncf=subset(cncf,chr<23)
  dipLogR <- fit$dipLogR
  
  #layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6), ncol=1))
  layout(matrix(c(1,1,2,2,3,3,4,4), ncol=1))
  par(mar=c(3,3,1,1), mgp=c(2, 0.7, 0))

  chr=mat$chrom
  len=table(chr)
  altcol=rep(c("light blue","gray"),12)[-c(23:24)]
  chr.col=rep(altcol,len)
  nmark=cncf$num.mark
  tmp=cumsum(len)
  start=c(1,tmp[-22]+1)
  end=tmp
  mid=start+len/2
  
  
  plot(mat$cnlr,pch=".",axes=F,cex=1.5,ylim=c(-3,3),col=c("grey","lightblue")[1+rep(cncf$chr-2*floor(cncf$chr/2),cncf$num.mark)],ylab="log-ratio",xlab=unique(orig$IGV$ID))
  abline(h=dipLogR, col="gray")
  points(rep(cncf$cnlr.median, cncf$num.mark), pch=".", cex=2, col="brown")
  axis(side=1,at=mid,1:22,cex.axis=1,las=2)
  axis(side=2)
  box()
  
  plot(mat$valor,axes=F,pch=".",cex=1.5,col=c("grey","lightblue")[1+rep(cncf$chr-2*floor(cncf$chr/2),cncf$num.mark)],ylab="log-odds-ratio",xlab="",ylim=c(-4,4))
  points(rep(sqrt(abs(cncf$mafR)), cncf$num.mark), pch=".", cex=2, col="brown")
  points(-rep(sqrt(abs(cncf$mafR)), cncf$num.mark), pch=".", cex=2, col="brown")
  axis(side=1,at=mid,1:22,cex.axis=1,las=2)
  axis(side=2)
  box()
   

	labeled=paste(
		"Purity = ",sprintf("%.3f",fit$purity),
		", Ploidy = ",sprintf("%.3f",fit$ploidy),
		", WGD = ",ifelse(fit$dipt>2,"Yes","No"),sep=""
	)
  plot(rep(cncf$cf.em,cncf$num.mark),axes=F,pch=".",cex=2,xlab="",ylab="Cellular fraction (EM)",ylim=c(0,1),main=labeled)
  axis(side=1,at=mid,1:22,cex.axis=1,las=2)
  axis(side=2,cex=0.8)
  box()
  abline(v=start,lty=3,col="gray")
  abline(v=end,lty=3,col="gray")
  #abline(h=c(0.2,0.4,0.6,0.8),lty=3,col="gray")
    
  
  # scale tcn so that very high copy numbers don't take up space
  tcnscaled <- cncf$tcn.em
  tcnscaled[cncf$tcn.em > 5 & !is.na(cncf$tcn.em)] = (5 + (tcnscaled[cncf$tcn.em > 5& !is.na(cncf$tcn.em)] -5)/3)
  matplot(cbind(rep(tcnscaled, cncf$num.mark), rep(cncf$lcn.em,cncf$num.mark)-0.1), pch=".", cex=3, col=1:2, lwd=1, ylab="Integer copy number (EM)", yaxt="n",xaxt="n")
  axis(2, at=c(0:5,5+(1:35)/3), labels=0:40)
  axis(side=1,at=mid,1:22,cex.axis=1,las=2)
  box()
  abline(v=start,lty=3,col="gray")
  abline(v=end,lty=3,col="gray")
  abline(h=c(0:5,5+(1:35)/3),lty=3,col="gray")
}
