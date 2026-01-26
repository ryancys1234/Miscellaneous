library(cvTools); library(glmnet); library(readr)
hpfs_met <- read_csv("~/165993/epi293/EPI293_TraitData/hpfs_met.csv")
hpfs_pheno <- read_csv("~/165993/epi293/EPI293_TraitData/hpfs_pheno.csv")
hpfs_pheno_covar <- read.delim("~/165993/epi293/EPI293_TraitData/hpfs_pheno_covar.txt")

hpfs <- merge(hpfs_met, hpfs_pheno, by.x = "id", by.y = "id")
hpfs <- hpfs[!is.na(hpfs$BMIcont),]
x <- names(hpfs)[2:289] # Metabolites
hpfs[, x] <- scale(hpfs[, x])
hpfs$T2D <- 0
hpfs[which(hpfs$diabetes == 1|hpfs$t2dbasediag == 1), "T2D"] <- 1
hpfs <- hpfs[, c(x, "ageyr", "BMIcont", "fast", "id", "smoking", "T2D")]
hpfs$predictScore <- 0

k <- 10
folds <- cvFolds(dim(hpfs)[1], K = k)
coeffMatrix <- matrix(0, nrow = length(x) + 1, ncol = k)
row.names(coeffMatrix) <- c("intercept", x)
for (i in 1:k) {
  fit <- cv.glmnet(as.matrix(hpfs[folds$subsets[folds$which != i], x]), hpfs[folds$subsets[folds$which != i], "BMIcont"], alpha = 0.5, family = "gaussian", type.measure = "mse")
  prediction <- predict(fit, as.matrix(hpfs[folds$subsets[folds$which == i], x]), s = "lambda.min", type = "response")
  coeffMatrix[, i] <- as.vector(coef(fit))
  hpfs[folds$subsets[folds$which == i], "predictScore"] = as.vector(prediction[, 1])
}
plot(hpfs[, "BMIcont"], hpfs[, "predictScore"])
cor(hpfs[, "predictScore"], hpfs[, "BMIcont"])
write.csv(hpfs, file = "~/EPI293/hpfs.csv", row.names = F) # Use for regression

colnames(hpfs)[which(colnames(hpfs) == "id")] <- "IID"
merged <- merge(x = hpfs_pheno_covar, y = hpfs[, c("IID", "predictScore")], by = "IID")
merged <- merged[, c(2, 1, 3:ncol(merged))]
write.table(merged, file = "~/EPI293/hpfs_pheno_covar.txt", sep = "\t", row.names = F, quote = F)

# coeffMatrix <- as.data.frame(coeffMatrix)
# coeffMatrix$AveCoeff <- 0
# coeffMatrix$Count <- 0
# for (i in 1:dim(coeffMatrix)[1]) {
#   coeffMatrix[i, "AveCoeff"] <- mean(t(coeffMatrix[i, 1:k]))
#   coeffMatrix[i, "Count"] <- sum(coeffMatrix[i, 1:k] != 0)
# }

# barplot(table(coeffMatrix$Count)); title(main = "Distribution of selected metabolites")
# hist(coeffMatrix$AveCoeff[coeffMatrix$AveCoeff < 20], breaks = 40, main = "Histogram of metabolite coefficients", xlab = "Metabolite coefficients")
# hist(coeffMatrix$AveCoeff[coeffMatrix$AveCoeff != 0 & coeffMatrix$AveCoeff < 20], breaks = 40, main = "Histogram of selected metabolite coefficients", xlab = "Selected metabolite coefficients")
# head(sort(coeffMatrix$AveCoeff, decreasing = T), n = 5)
# head(sort(coeffMatrix$AveCoeff), n = 5)
# colSums(coeffMatrix[, 1:10] != 0)
