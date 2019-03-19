plotSampleCNCF.custom<-function (jointseg, out, fit,main="")
{
    mat = jointseg
    mat = subset(mat, chrom < 23)
    out = subset(out, chrom < 23) #from chr to chrom
    cncf = fit$cncf
    cncf = subset(cncf, chrom < 23) #from chr to chrom
    dipLogR <- fit$dipLogR
    layout(matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6), ncol = 1))
    par(mar = c(3, 3, 1, 1), mgp = c(2, 0.7, 0))
    chr = mat$chrom
    len = table(chr)
    altcol = rep(c("light blue", "gray"), 12)[-c(23:24)]
    chr.col = rep(altcol, len)
    nmark = cncf$num.mark
    tmp = cumsum(len)
    start = c(1, tmp[-22] + 1)
    end = tmp
    mid = start + len/2
    plot(mat$cnlr, pch = ".", axes = F, cex = 1.5, ylim = c(-3,
        3), col = c("grey", "lightblue")[1 + rep(cncf$chr - 2 *
        floor(cncf$chr/2), cncf$num.mark)], ylab = "log-ratio", main=main)
    abline(h = dipLogR, col = "gray")
    points(rep(cncf$cnlr.median, cncf$num.mark), pch = ".", cex = 2,
        col = "brown")
    axis(side = 1, at = mid, 1:22, cex.axis = 1, las = 2)
    axis(side = 2)
    box()
    plot(mat$valor, axes = F, pch = ".", cex = 1.5, col = c("grey",
        "lightblue")[1 + rep(cncf$chr - 2 * floor(cncf$chr/2),
        cncf$num.mark)], ylab = "log-odds-ratio", ylim = c(-4,
        4))
    points(rep(sqrt(abs(cncf$mafR)), cncf$num.mark), pch = ".",
        cex = 2, col = "brown")
    points(-rep(sqrt(abs(cncf$mafR)), cncf$num.mark), pch = ".",
        cex = 2, col = "brown")
    axis(side = 1, at = mid, 1:22, cex.axis = 1, las = 2)
    axis(side = 2)
    box()
    plot(rep(cncf$cf.em, cncf$num.mark), axes = F, pch = ".",
        cex = 2, xlab = "Chromosome", ylab = "Cellular fraction (EM)",
        ylim = c(0, 1))
    axis(side = 1, at = mid, 1:22, cex.axis = 1, las = 2)
    axis(side = 2, cex = 0.8)
    box()
    abline(v = start, lty = 3, col = "gray")
    abline(v = end, lty = 3, col = "gray")
    tcnscaled <- cncf$tcn.em
    tcnscaled[cncf$tcn.em > 5 & !is.na(cncf$tcn.em)] = (5 + (tcnscaled[cncf$tcn.em >
        5 & !is.na(cncf$tcn.em)] - 5)/3)
    matplot(cbind(rep(tcnscaled, cncf$num.mark), rep(cncf$lcn.em,
        cncf$num.mark) - 0.1), pch = ".", cex = 3, col = 1:2,
        lwd = 1, ylab = "Integer copy number (EM)", yaxt = "n",
        xaxt = "n")
    axis(2, at = c(0:5, 5 + (1:35)/3), labels = 0:40)
    axis(side = 1, at = mid, 1:22, cex.axis = 1, las = 2)
    box()
    abline(v = start, lty = 3, col = "gray")
    abline(v = end, lty = 3, col = "gray")
    abline(h = c(0:5, 5 + (1:35)/3), lty = 3, col = "gray")
    plot(rep(cncf$cf, cncf$num.mark), axes = F, pch = ".", cex = 2,
        xlab = "Chromosome", ylab = "Cellular fraction (cncf)",
        ylim = c(0, 1))
    axis(side = 1, at = mid, 1:22, cex.axis = 1, las = 2)
    axis(side = 2, cex = 0.8)
    box()
    abline(v = start, lty = 3, col = "gray")
    abline(v = end, lty = 3, col = "gray")
    tcnscaled <- cncf$tcn
    tcnscaled[cncf$tcn > 5 & !is.na(cncf$tcn)] = (5 + (tcnscaled[cncf$tcn >
        5 & !is.na(cncf$tcn)] - 5)/3)
    matplot(cbind(rep(tcnscaled, cncf$num.mark), rep(cncf$lcn,
        cncf$num.mark) - 0.1), pch = ".", cex = 3, col = 1:2,
        lwd = 1, ylab = "Integer copy number (cncf)", yaxt = "n",
        xaxt = "n")
    axis(2, at = c(0:5, 5 + (1:35)/3), labels = 0:40)
    axis(side = 1, at = mid, 1:22, cex.axis = 1, las = 2)
    box()
    abline(v = start, lty = 3, col = "gray")
    abline(v = end, lty = 3, col = "gray")
    abline(h = c(0:5, 5 + (1:35)/3), lty = 3, col = "gray")
}
