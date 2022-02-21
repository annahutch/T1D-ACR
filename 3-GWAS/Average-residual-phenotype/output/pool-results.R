## code to merge results
library(data.table)

imp1 <- fread("imp1-out.assoc.txt")
imp2 <- fread("imp2-out.assoc.txt")
imp3 <- fread("imp3-out.assoc.txt")
imp4 <- fread("imp4-out.assoc.txt")
imp5 <- fread("imp5-out.assoc.txt")

allimpres <- list(imp1, imp2, imp3, imp4, imp5)

#newdf <- NULL

#for(i in 1:7819188){
  
#  newdf[[i]] <- data.frame(rbindlist(lapply(allimpres,function(x) x[i,9])), rbindlist(lapply(allimpres,function(x) x[i,10])), rbindlist(lapply(allimpres,function(x) x[i,11])))
  
#}

imp1.split <- split(allimpres[[1]][,9:11], seq(nrow(allimpres[[1]][,9:11])))
imp2.split <- split(allimpres[[2]][,9:11], seq(nrow(allimpres[[2]][,9:11])))
imp3.split <- split(allimpres[[3]][,9:11], seq(nrow(allimpres[[3]][,9:11])))
imp4.split <- split(allimpres[[4]][,9:11], seq(nrow(allimpres[[4]][,9:11])))
imp5.split <- split(allimpres[[5]][,9:11], seq(nrow(allimpres[[5]][,9:11])))

newdf <- mapply(rbind, imp1.split, imp2.split, imp3.split, imp5.split, imp5.split, SIMPLIFY=FALSE)

saveRDS(newdf, "seperate-imp-lists.RDS")

# https://bookdown.org/mwheymans/bookmi/rubins-rules.html
# https://stats.stackexchange.com/questions/327237/calculating-pooled-p-values-manually

out <- lapply(newdf, function(x){
  
  mat <- x
  mat <- as.data.frame(mat)
  
  (pooledMean <- mean(mat[,1]))
  
  (betweenVar <- mean(mat[,2]^2)) # mean of variances
  (withinVar <- sd(mat[,1])^2) # variance of variances
  (dfCorrection <- (nrow(mat)+1)/(nrow(mat))) # dfCorrection
  
  (totVar <- betweenVar + withinVar*dfCorrection) # total variance
  (pooledSE <- sqrt(totVar)) # standard error
  
  lambda <- (withinVar + (withinVar/5))/totVar # fraction of missing information
  n <- 1739
  k <- 2
  nu_old <- (5-1)/lambda^2  
  nu_com <- n-k
  nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
  (nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
  
  pooledP <- pt(q = abs(pooledMean / pooledSE), df = nu_BR, lower.tail = FALSE) * 2
  
  data.frame(pooledBeta = pooledMean, pooledSE, pooledP)
  
})

final.tmp <- rbindlist(out)
final.tmp$chr <- imp1$chr
final.tmp$ps <- imp1$ps
final.tmp$rs <- imp1$rs
final.tmp$allele1 <- imp1$allele1
final.tmp$allele0 <- imp1$allele0
final.tmp$af <- imp1$af

saveRDS(final.tmp, "pooled-res.RDS")

