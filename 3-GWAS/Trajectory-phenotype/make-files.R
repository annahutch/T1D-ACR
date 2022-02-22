merged_fam <- read.table("../../1-GenotypeQCImp/Merged/merged.fam", header = F)
matchIDs <- readRDS("../../1-GenotypeQCImp/sample-keys/genotypedsamples-IDs.RDS")

phenotypes.imp1 <- readRDS("../../2-PhenotypeDataPrep/c-MakePhenotypes/Trajectory-pheno/phenotypes-iter1.RDS")
phenotypes.imp2 <- readRDS("../../2-PhenotypeDataPrep/c-MakePhenotypes/Trajectory-pheno/phenotypes-iter2.RDS")
phenotypes.imp3 <- readRDS("../../2-PhenotypeDataPrep/c-MakePhenotypes/Trajectory-pheno/phenotypes-iter3.RDS")
phenotypes.imp4 <- readRDS("../../2-PhenotypeDataPrep/c-MakePhenotypes/Trajectory-pheno/phenotypes-iter4.RDS")
phenotypes.imp5 <- readRDS("../../2-PhenotypeDataPrep/c-MakePhenotypes/Trajectory-pheno/phenotypes-iter5.RDS")

merged_fam$altID <- matchIDs$pheno_ID[match(merged_fam$V2, matchIDs$geno_ID)]

# imp 1
merged_fam$V6 <- phenotypes.imp1$brms.Slope[match(merged_fam$altID, phenotypes.imp1$StudyID)]
write.table(merged_fam[,-7], file = "imp1/merged.fam", row.names = F, col.names = F, quote = F)
system("cp ../../1-GenotypeQCImp/Merged/merged.bim imp1/") # copy over bim and bed files
system("cp ../../1-GenotypeQCImp/Merged/merged.bed imp1/") 

# imp 2
merged_fam$V6 <- phenotypes.imp2$brms.Slope[match(merged_fam$altID, phenotypes.imp2$StudyID)]
write.table(merged_fam[,-7], file = "imp2/merged.fam", row.names = F, col.names = F, quote = F)
system("cp ../../1-GenotypeQCImp/Merged/merged.bim imp2/") # copy over bim and bed files
system("cp ../../1-GenotypeQCImp/Merged/merged.bed imp2/") 

# imp 3
merged_fam$V6 <- phenotypes.imp3$brms.Slope[match(merged_fam$altID, phenotypes.imp3$StudyID)]
write.table(merged_fam[,-7], file = "imp3/merged.fam", row.names = F, col.names = F, quote = F)
system("cp ../../1-GenotypeQCImp/Merged/merged.bim imp3/") # copy over bim and bed files
system("cp ../../1-GenotypeQCImp/Merged/merged.bed imp3/") 

# imp 4
merged_fam$V6 <- phenotypes.imp4$brms.Slope[match(merged_fam$altID, phenotypes.imp4$StudyID)]
write.table(merged_fam[,-7], file = "imp4/merged.fam", row.names = F, col.names = F, quote = F)
system("cp ../../1-GenotypeQCImp/Merged/merged.bim imp4/") # copy over bim and bed files
system("cp ../../1-GenotypeQCImp/Merged/merged.bed imp4/") 

# imp 5
merged_fam$V6 <- phenotypes.imp5$brms.Slope[match(merged_fam$altID, phenotypes.imp5$StudyID)]
write.table(merged_fam[,-7], file = "imp5/merged.fam", row.names = F, col.names = F, quote = F)
system("cp ../../1-GenotypeQCImp/Merged/merged.bim imp5/") # copy over bim and bed files
system("cp ../../1-GenotypeQCImp/Merged/merged.bed imp5/") 
