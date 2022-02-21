library(GenomicRanges)
library(karyoploteR)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(corrcoverage)

traj_res <- readRDS("../../3-GWAS/Trajectory-phenotype/output/combine-results/final-pooled-res.RDS")

chr_range <- 5
bp_range <- c(20409042 - 200000, 20409042 + 200000)

traj_res_myregion_tmp <- traj_res[which(traj_res$chr == chr_range),]
traj_res_myregion <- traj_res_myregion_tmp[which(traj_res_myregion_tmp$ps >= bp_range[1] & traj_res_myregion_tmp$ps <= bp_range[2]),]
traj_res_myregion[ , chr := paste0('chr', chr)]

Z <- traj_res_myregion$pooledBeta/traj_res_myregion$pooledSE

traj_res_myregion$PP <- ppfunc(z = Z, V = traj_res_myregion$pooledSE^2)

colnames(traj_res_myregion)[c(1,2,9)] <- c("CHR19","BP19","P")
trajGRanges <- makeGRangesFromDataFrame(data.frame(traj_res_myregion), start.field = 'BP19', end.field = 'BP19', seqnames.field = 'CHR19', ignore.strand = T, keep.extra.columns = T)

png("CDH18.png", width = 600, height = 500)

cex.label <- 1
ymax <- 10
label.margin <- 0.05
# The % of the total plotting space to leave in each margin between the tracks
autotrack.margin <- 0.18

pp <- getDefaultPlotParams(plot.type=4)

pp$leftmargin <- 0.1

kp <- plotKaryotype(plot.type = 4, labels.plotter = NULL, zoom = paste0("chr",chr_range,":",bp_range[1],"-",bp_range[2]), plot.params = pp)
kpAddBaseNumbers(kp, add.units = T, cex = 0.75, tick.dist = 1e6)

title(main = 'CDH18 region', cex.main = 1.5)

at<-autotrack(3, 3, margin = autotrack.margin)
kp <- kpPlotManhattan(kp, data = trajGRanges, pval = trajGRanges$PP, points.col  =  "black", r0 = at$r0, r1 = at$r1, ymax = 1, points.cex=1, logp = FALSE, suggestive.col = "white", genomewide.col = "white")
kpAddLabels(kp, labels  =  "PP", srt = 90, pos = 3, r0 = at$r0, r1 = at$r1, cex = cex.label, label.margin = 0.06)
kpAxis(kp, ymin = 0, ymax = 1, r0 = at$r0, r1 = at$r1)

kpPoints(kp, data=trajGRanges[which.min(trajGRanges$P),], y = trajGRanges$PP[which.min(trajGRanges$P)], pch=1, cex=2, col="red", ymax=1, r0=at$r0, r1 = at$r1)
kpText(kp, data=trajGRanges[which.min(trajGRanges$P),], y = trajGRanges$PP[which.min(trajGRanges$P)], labels = "rs145715205", pos=2, ymax=1, r0=at$r0, col = "red", r1 = at$r1)

at<-autotrack(2, 3, margin = 0.07)
kp <- kpPlotManhattan(kp, data = trajGRanges, pval = trajGRanges$P, points.col  =  "black", r0 = at$r0, r1 = at$r1, ymax = ymax, points.cex=1)
kpAddLabels(kp, labels  =  expression(-log[10] ~ (p)), srt = 90, pos = 3, r0 = at$r0, r1 = at$r1, cex = cex.label, label.margin = 0.06)
kpAxis(kp, ymin = 0, ymax = ymax, r0 = at$r0, r1 = at$r1)

kpPoints(kp, data=trajGRanges[which.min(trajGRanges$P),], y = -log10(trajGRanges$P[which.min(trajGRanges$P)]), pch=1, cex=2, col="red", ymax=ymax, r0=at$r0, r1 = at$r1)

kpText(kp, data=trajGRanges[which.min(trajGRanges$P),], y = -log10(trajGRanges$P[which.min(trajGRanges$P)]), labels = "rs145715205", pos=2, ymax=ymax, r0=at$r0, col = "red", r1 = 0.625)

# Needed to install org.Hs.eg.db with BiocManager::install to get this running
genes.data <- makeGenesDataFromTxDb(karyoplot = kp, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene)

genes.data <- addGeneNames(genes.data)
genes.data.merged <- mergeTranscripts(genes.data)
kp <- kpPlotGenes(kp, data = genes.data.merged, r0 = 0, r1 = 0.2, gene.name.cex = 0.8, gene.name.position = 'right')

dev.off()
