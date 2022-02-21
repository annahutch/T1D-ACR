merged_fam <- read.table("../../1-GenotypeQCImp/Merged/merged.fam", header = F)
matchIDs <- readRDS("../../1-GenotypeQCImp/sample-keys/genotypedsamples-IDs.RDS")

phenotypes <- readRDS("../../2-PhenotypeDataPrep/c-MakePhenotypes/Secondary-phenotypes/secondary-phenotypes.RDS ")

merged_fam$altID <- matchIDs$pheno_ID[match(merged_fam$V2, matchIDs$geno_ID)]

# imp 1
merged_fam$V6 <- phenotypes$latestACR.imp1[match(merged_fam$altID, phenotypes$StudyID)]
write.table(merged_fam[,-7], file = "imp1/merged.fam", row.names = F, col.names = F, quote = F)
system("cp ../../1-GenotypeQCImp/Merged/merged.bim imp1/") # copy over bim and bed files
system("cp ../../1-GenotypeQCImp/Merged/merged.bed imp1/") 

# imp 2
merged_fam$V6 <- phenotypes$latestACR.imp2[match(merged_fam$altID, phenotypes$StudyID)]
write.table(merged_fam[,-7], file = "imp2/merged.fam", row.names = F, col.names = F, quote = F)
system("cp ../../1-GenotypeQCImp/Merged/merged.bim imp2/") # copy over bim and bed files
system("cp ../../1-GenotypeQCImp/Merged/merged.bed imp2/") 

# imp 3
merged_fam$V6 <- phenotypes$latestACR.imp3[match(merged_fam$altID, phenotypes$StudyID)]
write.table(merged_fam[,-7], file = "imp3/merged.fam", row.names = F, col.names = F, quote = F)
system("cp ../../1-GenotypeQCImp/Merged/merged.bim imp3/") # copy over bim and bed files
system("cp ../../1-GenotypeQCImp/Merged/merged.bed imp3/") 

# imp 4
merged_fam$V6 <- phenotypes$latestACR.imp4[match(merged_fam$altID, phenotypes$StudyID)]
write.table(merged_fam[,-7], file = "imp4/merged.fam", row.names = F, col.names = F, quote = F)
system("cp ../../1-GenotypeQCImp/Merged/merged.bim imp4/") # copy over bim and bed files
system("cp ../../1-GenotypeQCImp/Merged/merged.bed imp4/") 

# imp 5
merged_fam$V6 <- phenotypes$latestACR.imp5[match(merged_fam$altID, phenotypes$StudyID)]
write.table(merged_fam[,-7], file = "imp5/merged.fam", row.names = F, col.names = F, quote = F)
system("cp ../../1-GenotypeQCImp/Merged/merged.bim imp5/") # copy over bim and bed files
system("cp ../../1-GenotypeQCImp/Merged/merged.bed imp5/") 
