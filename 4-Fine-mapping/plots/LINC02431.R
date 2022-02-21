library(GenomicRanges)
library(karyoploteR)
library(data.table)
library(corrcoverage)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

latestACR_res <- readRDS("../../3-GWAS/Latest-residual-phenotype/output/combine-results/final-pooled-res.RDS")

chr_range <- 4
bp_range <- c(171934268 - 200000, 171934268 + 200000)

latestACR_res_myregion_tmp <- latestACR_res[which(latestACR_res$chr == chr_range),]
latestACR_res_myregion <- latestACR_res_myregion_tmp[which(latestACR_res_myregion_tmp$ps >= bp_range[1] & latestACR_res_myregion_tmp$ps <= bp_range[2]),]
latestACR_res_myregion[ , chr := paste0('chr', chr)]

Z <- latestACR_res_myregion$pooledBeta/latestACR_res_myregion$pooledSE

latestACR_res_myregion$PP <- ppfunc(z = Z, V = latestACR_res_myregion$pooledSE^2)

credset(latestACR_res_myregion$PP, thr = 0.95)

colnames(latestACR_res_myregion)[c(1,2,9)] <- c("CHR19","BP19","P")
latestACRGRanges <- makeGRangesFromDataFrame(data.frame(latestACR_res_myregion), start.field = 'BP19', end.field = 'BP19', seqnames.field = 'CHR19', ignore.strand = T, keep.extra.columns = T)


png("LINC02431.png", width = 600, height = 500)

cex.label <- 1
ymax <- 10
label.margin <- 0.05
# The % of the total plotting space to leave in each margin between the tracks
autotrack.margin <- 0.18

pp <- getDefaultPlotParams(plot.type=4)

pp$leftmargin <- 0.1

kp <- plotKaryotype(plot.type = 4, labels.plotter = NULL, zoom = paste0("chr",chr_range,":",bp_range[1],"-",bp_range[2]), plot.params = pp)
kpAddBaseNumbers(kp, add.units = T, cex = 0.75, tick.dist = 1e6)

title(main = 'LINC02431 region', cex.main = 1.5)

at<-autotrack(3, 3, margin = autotrack.margin)
kp <- kpPlotManhattan(kp, data = latestACRGRanges, pval = latestACRGRanges$PP, points.col  =  "black", r0 = at$r0, r1 = at$r1, ymax = 1, points.cex=1, logp = FALSE, suggestive.col = "white", genomewide.col = "white")
kpAddLabels(kp, labels  =  "PP", srt = 90, pos = 3, r0 = at$r0, r1 = at$r1, cex = cex.label, label.margin = 0.06)
kpAxis(kp, ymin = 0, ymax = 1, r0 = at$r0, r1 = at$r1)

kpPoints(kp, data=latestACRGRanges[which.min(latestACRGRanges$P),], y = latestACRGRanges$PP[which.min(latestACRGRanges$P)], pch=1, cex=2, col="red", ymax=1, r0=at$r0, r1 = at$r1)
kpText(kp, data=latestACRGRanges[which.min(latestACRGRanges$P),], y = latestACRGRanges$PP[which.min(latestACRGRanges$P)], labels = "rs150766792", pos=2, ymax=1, r0=at$r0, col = "red", r1 = at$r1)

# add blue circles around other SNPs in credset
kpPoints(kp, data=latestACRGRanges[which(latestACRGRanges$rs=="4:171864867:G:A"),], y = latestACRGRanges$PP[which(latestACRGRanges$rs=="4:171864867:G:A")], pch=1, cex=2, col="blue", ymax=1, r0=at$r0, r1 = at$r1)
kpPoints(kp, data=latestACRGRanges[which(latestACRGRanges$rs=="4:171957178:C:G"),], y = latestACRGRanges$PP[which(latestACRGRanges$rs=="4:171957178:C:G")], pch=1, cex=2, col="blue", ymax=1, r0=at$r0, r1 = at$r1)

at<-autotrack(2, 3, margin = 0.07)
kp <- kpPlotManhattan(kp, data = latestACRGRanges, pval = latestACRGRanges$P, points.col  =  "black", r0 = at$r0, r1 = at$r1, ymax = ymax, points.cex=1)
kpAddLabels(kp, labels  =  expression(-log[10] ~ (p)), srt = 90, pos = 3, r0 = at$r0, r1 = at$r1, cex = cex.label, label.margin = 0.06)
kpAxis(kp, ymin = 0, ymax = ymax, r0 = at$r0, r1 = at$r1)

kpPoints(kp, data=latestACRGRanges[which.min(latestACRGRanges$P),], y = -log10(latestACRGRanges$P[which.min(latestACRGRanges$P)]), pch=1, cex=2, col="red", ymax=ymax, r0=at$r0, r1 = at$r1)

kpText(kp, data=latestACRGRanges[which.min(latestACRGRanges$P),], y = -log10(latestACRGRanges$P[which.min(latestACRGRanges$P)]), labels = "rs150766792", pos=2, ymax=ymax, r0=at$r0, col = "red", r1 = at$r1)

kpPoints(kp, data=latestACRGRanges[which(latestACRGRanges$rs=="4:171864867:G:A"),], y = -log10(latestACRGRanges$P[which(latestACRGRanges$rs=="4:171864867:G:A")]), pch=1, cex=2, col="blue", ymax=ymax, r0=at$r0, r1 = at$r1)
kpPoints(kp, data=latestACRGRanges[which(latestACRGRanges$rs=="4:171957178:C:G"),], y = -log10(latestACRGRanges$P[which(latestACRGRanges$rs=="4:171957178:C:G")]), pch=1, cex=2, col="blue", ymax=ymax, r0=at$r0, r1 = at$r1)

# Needed to install org.Hs.eg.db with BiocManager::install to get this running
genes.data <- makeGenesDataFromTxDb(karyoplot = kp, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene)

genes.data <- addGeneNames(genes.data)
genes.data.merged <- mergeTranscripts(genes.data)
kp <- kpPlotGenes(kp, data = genes.data.merged, r0 = 0, r1 = 0.2, gene.name.cex = 0.8, gene.name.position = 'right')

dev.off()
