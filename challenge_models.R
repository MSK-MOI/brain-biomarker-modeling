library(bbutilities)
options(digits=2)


# Build models for all 3 challenges and report statistics
filenames1 <- c("sc1_Phase1_GE_Phenotype.tsv", "sc1_Phase1_GE_Outcome.tsv", "sc1_Phase1_GE_FeatureMatrix.tsv")
ld1 <- pipeline(filenames1, threshold=12, cores=6, write="subchallenge1")
filenames2 <- c("sc2_Phase1_CN_Phenotype.tsv", "sc2_Phase1_CN_Outcome.tsv", "sc2_Phase1_CN_FeatureMatrix.tsv")
ld2 <- pipeline(filenames2, threshold=8, cores=6, write="subchallenge2")
filenames3 <- c("sc3_Phase1_CN_GE_Phenotype.tsv", "sc3_Phase1_CN_GE_Outcome.tsv", "sc3_Phase1_CN_GE_FeatureMatrix.tsv")
ld3 <- pipeline(filenames3, threshold=7, cores=6, write="subchallenge3")


# Use the following (after the above) to print the reports quickly for selected values of the threshold (does not repeat the main calculations)
ld1 <- pipeline(ld=ld1, merged_ranking=ld1[[3]], filenames1, threshold=12, cores=6, write="subchallenge1")
ld2 <- pipeline(ld=ld2, merged_ranking=ld2[[3]], filenames2, threshold=8, cores=6, write="subchallenge2")
ld3 <- pipeline(ld=ld3, merged_ranking=ld3[[3]], filenames3, threshold=7, cores=6, write="subchallenge3")


# Use the following to calculate a summary statistic of the cross-validation efficacy for multiple thresholds on the number of features.
stats1 <- lapply(1:40, function(j) {
        ld1 <- pipeline(ld=ld1, merged_ranking=ld1[[3]], filenames1, threshold=j, cores=6, write="subchallenge1_experiments")
        return(ld1[[4]])
    })
mean_cv_aucs1 <- sapply(1:40, function(j){return(stats1[[j]][1])})
accs1 <- sapply(1:40, function(j){return(stats1[[j]][2])})
specs1 <- sapply(1:40, function(j){return(stats1[[j]][4])})
cvaccs1 <- sapply(1:40, function(j){return(stats1[[j]][5])})

stats2 <- lapply(1:40, function(j) {
        ld2 <- pipeline(ld=ld2, merged_ranking=ld2[[3]], filenames2, threshold=j, cores=6, write="subchallenge2_experiments")
        return(ld2[[4]])
    })
mean_cv_aucs2 <- sapply(1:40, function(j){return(stats2[[j]][1])})
accs2 <- sapply(1:40, function(j){return(stats2[[j]][2])})
specs2 <- sapply(1:40, function(j){return(stats2[[j]][4])})
cvaccs2 <- sapply(1:40, function(j){return(stats2[[j]][5])})

stats3 <- lapply(1:40, function(j) {
        ld3 <- pipeline(ld=ld3, merged_ranking=ld3[[3]], filenames3, threshold=j, cores=6, write="subchallenge3_experiments")
        return(ld3[[4]])
    })
mean_cv_aucs3 <- sapply(1:40, function(j){return(stats3[[j]][1])})
accs3 <- sapply(1:40, function(j){return(stats3[[j]][2])})
specs3 <- sapply(1:40, function(j){return(stats3[[j]][4])})
cvaccs3 <- sapply(1:40, function(j){return(stats3[[j]][5])})


# Accuracy and specificity versus CV AUC
par(mfrow=c(1,3))
plot(c(1:40), mean_cv_aucs1, pch=16, xlim=c(0,40), ylim=c(0,1), xlab="threshold, challenge 1", ylab="CV AUCs (black) and final accuracy (blue) and final specificity (green)")
points(cbind(1:40, accs1), pch=16, col="blue")
points(cbind(1:40, specs1), pch=16, col="green")
points(cbind(1:40, cvaccs1), pch=16, col="red")

plot(c(1:40), mean_cv_aucs2, pch=16, xlim=c(0,40), ylim=c(0,1), xlab="threshold, challenge 2", ylab="CV AUCs (black) and final accuracy (blue) and final specificity (green)")
points(cbind(1:40, accs2), pch=16, col="blue")
points(cbind(1:40, specs2), pch=16, col="green")
points(cbind(1:40, cvaccs2), pch=16, col="red")

plot(c(1:40), mean_cv_aucs3, pch=16, xlim=c(0,40), ylim=c(0,1), xlab="threshold, challenge 3", ylab="CV AUCs (black) and final accuracy (blue) and final specificity (green)")
points(cbind(1:40, accs3), pch=16, col="blue")
points(cbind(1:40, specs3), pch=16, col="green")
points(cbind(1:40, cvaccs3), pch=16, col="red")


