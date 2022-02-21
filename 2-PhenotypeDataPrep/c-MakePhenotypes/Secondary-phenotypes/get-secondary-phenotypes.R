library(mice) 
library(dplyr)

impres <- readRDS("../../b-MultipleImputations/multipleimpres.RDS")
forstudyid <- readRDS("../../a-CleanData/final_df.RDS")
forstudyid <- forstudyid[-which(is.na(forstudyid$treatment.allocation)),]

phenodf <- data.frame(StudyID = forstudyid$StudyID, latestACR.imp1 = NA, latestACR.imp2 = NA, latestACR.imp3 = NA, latestACR.imp4 = NA, latestACR.imp5 = NA, averes.imp1 = NA, averes.imp2 = NA, averes.imp3 = NA, averes.imp4 = NA, averes.imp5 = NA)

# for each individual, extract the residual from the most recent visit and the average residual value across all visits
for(i in 1:5){
  
  df <- complete(impres, i)
  
  df <- df[-which(is.na(df$treatment.allocation)),]
  
  df$sex <- as.factor(df$sex)
  df$cohort <- as.factor(df$cohort)
  df$treatment.allocation <- as.factor(df$treatment.allocation)
  
  full_model <- lm(log(meanACR) ~ sex + visit_duration + visit_age + HbA1c + treatment.allocation + cohort + BMI + SBP + DBP, data = df)
  
  df$residual <- full_model$residuals
  
  # residual value at latest visit
  out <- df %>% group_by(StudyID) %>% mutate(latestACR = residual[which(max(visit_age, na.rm = TRUE)  == visit_age)[1]])
  phenodf[, i +1] <- out$latestACR
  
  # average residual
  averes <- aggregate(residual~StudyID, df, FUN=mean) 
  phenodf[, i + 6] <- averes$residual[match(df$StudyID, averes$StudyID)]
  
}

# keep only one row for each individual
new_phenodf <- phenodf[!duplicated(phenodf$StudyID),]

saveRDS(new_phenodf, "secondary-phenotypes.RDS")
