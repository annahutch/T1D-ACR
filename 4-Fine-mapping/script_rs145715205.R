library(susieR)
library(data.table)
library(tidyr)
library(rstatix)
library(plyr)
library(dplyr)
library(biomaRt)
library(corrcoverage)

args <- "5:20409042:A:C"

GWAS_data <- readRDS("../3-GWAS/Trajectory-phenotype/output/combine-results/final-pooled-res.RDS")

index <- which(GWAS_data$rs==args)
indexSNP_chr <- GWAS_data$chr[index]
indexSNP_location <- GWAS_data$ps[index]

subset_chr <- GWAS_data[which(GWAS_data$chr==indexSNP_chr), ]

snpsinregions <- subset_chr[ subset_chr$ps >= indexSNP_location-200000 & subset_chr$ps <= indexSNP_location+200000, ]

write.table(snpsinregions$rs, file = paste0(args,"_region.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# need to update .bim files with chr:pos:a1:a2 naming convention to match snps
bimfile <- fread(paste0("LDdata/EUR_phase3.MAF_001.chr",indexSNP_chr,".bim"))
bimfile$V2 <- paste0(bimfile$V1,":",bimfile$V4,":",bimfile$V6,":",bimfile$V5)
write.table(bimfile, file = paste0("LDdata/EUR_phase3.MAF_001.chr",indexSNP_chr,".bim"), quote = FALSE, row.names = FALSE, col.names = FALSE)

system("plink --bfile LDdata/EUR_phase3.MAF_001.chr5 --extract 5:20409042:A:C_region.txt --r2 --ld-window-r2 0 --ld-window-kb 1000 --ld-window 9999 --out plink_output/5:20409042:A:C_ld")
#system("plink --bfile LDdata/EUR_phase3.MAF_001.chr22 --extract 22:50050517:C:G_region.txt --r2 --ld-window-r2 0 --ld-window-kb 1000 --ld-window 9999 --out plink_output/22:50050517:C:G_ld")

### make LD matrix
cor_ld_long <- fread(paste0("plink_output/", args, "_ld.ld"))
common_snps <- intersect(cor_ld_long$SNP_A, cor_ld_long$SNP_B)
cor_ld_long <- cor_ld_long %>% filter(SNP_A %in% common_snps)
cor_ld_long <- cor_ld_long %>% filter(SNP_B %in% common_snps)

var_file <- fread(paste0(args,"_region.txt"), header = F)

var_file <- var_file %>% filter(V1 %in% common_snps)
var_file$V2 <- var_file$V1
var_file$R <- 1

cor_ld_long$R <- sqrt(cor_ld_long$R2)

colnames(cor_ld_long)[c(3,6)] <- c("V1","V2")
cor_ld_long <- cor_ld_long[,c(3,6,8)]

cor_ld_long_2 <- rbind(var_file,cor_ld_long)

out <- pivot_wider(cor_ld_long_2,names_from = "V2" , values_from = "R", names_sort = F)

# https://community.rstudio.com/t/convert-tibble-into-matrix/68262
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}

my_matrix <- make_matrix(select(out,-V1), out$V1)

my_matrix[lower.tri(my_matrix)] <- t(my_matrix)[lower.tri(my_matrix)]

# ony include SNPs present in LD matrix
snpsinregions <- snpsinregions[-which(!snpsinregions$rs %in% colnames(my_matrix)),]

# convert p-vals to Z scores
zp = -qnorm(snpsinregions$pooledP/2)

fitted_rss <- susie_rss(z = zp, R = my_matrix, L = 10, coverage = 0.9)

saveRDS(fitted_rss, file = paste0("SuSiE-output/SuSiE-output-", args, ".RDS"))


##########

snpsinregions$Z <- snpsinregions$pooledBeta/snpsinregions$pooledSE

pps <- ppfunc(z = snpsinregions$Z, V = snpsinregions$pooledSE^2)

snpsinregions$PP <- pps

one <- plot(snpsinregions$ps, -log10(snpsinregions$pooledP))
two <- plot(snpsinregions$ps, snpsinregions$PP)

cs <- credset(pp = pps, thr = 0.95)
