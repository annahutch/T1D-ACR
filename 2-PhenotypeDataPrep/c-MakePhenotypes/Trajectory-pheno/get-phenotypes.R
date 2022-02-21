## to be run as an array job on HPC for each iteration of imputation

library(brms)
library(mice)

task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)
i <- as.numeric(task_id_string)

impres <- readRDS("../../b-MultipleImputations/multipleimpres.RDS")

df <- complete(impres, i)
# remove those with missing treatment allocation
df <- df[-which(is.na(df$treatment.allocation)),]

df$sex <- as.factor(df$sex)
df$cohort <- as.factor(df$cohort)
df$treatment.allocation <- as.factor(df$treatment.allocation)

full_model <- lm(log(meanACR) ~ sex + visit_duration + visit_age + HbA1c + treatment.allocation + cohort + BMI + SBP + DBP, data = df)

saveRDS(full_model, paste0("lm-model-iter",i,".RDS"))

# extract residuals
df$residual <- full_model$residuals

saveRDS(df, paste0("df-withresidual-iter",i,".RDS"))

# average residual = adjusted albumin excretion phenotype
averes <- aggregate(residual~StudyID, df, FUN=mean) 

# consider other summary measures - intercept and slope from linear model fitted to trajectory
lmres <- data.frame(StudyID = unique(df$StudyID), Intercept = rep(NA, times = length(unique(df$StudyID))), Slope = rep(NA, times = length(unique(df$StudyID))), novisits = rep(NA, times = length(unique(df$StudyID))))

for(i in 1:length(unique(df$StudyID))){
  lmres[i,2:3] = coef(lm(residual~visit_duration, data = df[which(df$StudyID==unique(df$StudyID)[i]),]))
  lmres[i,4] = length(which(df$StudyID==unique(df$StudyID)[i]))
}

# add column for average residual
lmres$AveResidual <- averes$residual[match(unique(lmres$StudyID), averes$StudyID)]

# model: variable intercept + variable slope
brms.model <- brm(residual ~ 1 + visit_duration + (1 + visit_duration | StudyID), iter = 8000, data = df)
saveRDS(brms.model, paste0("brmsmodel-iter",task_id,".RDS"))

# extract posterior samples
posteriorsamples <- posterior_samples(brms.model)

# extract mean of posterior samples for fixed intercept and fixed slope
fixedintercept <- mean(posteriorsamples$b_Intercept)
fixedslope <- mean(posteriorsamples$b_visit_duration)

# extract variable/random intercepts and slopes
intercept_samples <- posteriorsamples[grepl("Intercept", names(posteriorsamples), fixed = TRUE)][-c(1,2,3)] 
StudyID<- sub(",.*", "", sub(".*\\[", "", names(intercept_samples)))
lmres$brms.Intercept <- apply(intercept_samples, 2, mean)[match(lmres$StudyID, StudyID)] + fixedintercept

slope_samples <- posteriorsamples[grepl("visit_duration", names(posteriorsamples), fixed = TRUE)][-c(1,2,3)] 
StudyID <- sub(",.*", "", sub(".*\\[", "", names(slope_samples)))
lmres$brms.Slope <- apply(slope_samples, 2, mean)[match(lmres$StudyID, StudyID)] + fixedslope

####### Add Study IDs from phenotype and genotype data to help with matching up plink files
forstudyid <- readRDS("final_df.RDS")
forstudyid <- forstudyid[-which(is.na(forstudyid$treatment.allocation)),]

lmres$StudyID <- unique(forstudyid$StudyID)

saveRDS(lmres, paste0("phenotypes-iter",task_id,".RDS"))
