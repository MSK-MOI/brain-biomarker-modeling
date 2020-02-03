
#' @importFrom stats glm
#' @importFrom stats binomial
#' @importFrom stats predict
#' @importFrom stats t.test
#' @importFrom stats aov
#' @importFrom stats summary.aov
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom utils read.csv
#' @importFrom utils write.csv
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom methods slot
#' @import pbmcapply
#' @import parallel
#' @import ROCR
#' @import nnet
#' @import caret

library(stats)
library(utils)
library(pbmcapply)
library(parallel)
library(ROCR)
library(nnet)
library(caret)

#' Load data (from brain biomarker challenge files)
#'
#' @param filenames Vector of filenames: phenotype, outcome, and feature data.
#' @param drop.zeros Whether to drop features that are identically zero.
#' @return list(df, feature_names). df$type, df$outcome, df$gene1, df$gene2, etc.
#' @export
load_data <- function(filenames, drop.zeros=FALSE) {
    categories <- read.csv(filenames[1], header=T, stringsAsFactors=F, sep="\t", row.names=1)
    type <- categories[,c("WHO_GRADING","CANCER_TYPE"),drop=F]

    outcome <- read.table(filenames[2], header=T, stringsAsFactors=F)
    rownames(outcome) <- outcome$PATIENTID
    outcome <- outcome[,-1, drop=F]

    gv <- read.table(filenames[3], header=T, stringsAsFactors=F)
    rownames(gv) <- gv$PATIENTID
    gv <- gv[,-1]

    if(drop.zeros) {
        gv <- gv[,sapply(1:(dim(gv)[2]), FUN=function(j){!all(gv[,j]==0)})]
    }

    feature_names <- colnames(gv)

    if(!( all(rownames(type)==rownames(outcome)) && all(rownames(outcome)==rownames(gv)) )) {
        stop("Patient IDs don't all match.")
    }

    df <- data.frame(cbind(type$CANCER_TYPE, outcome$SURVIVAL_STATUS, gv, stringsAsFactors=F))
    colnames(df)[1] <- "CANCER_TYPE"
    colnames(df)[2] <- "SURVIVAL_STATUS"
    return(list(df, feature_names))
}


#' Run t-test to assess relevance of a variable
#' 
#' @param X Data frame with numerical data in first column, and SURVIVAL_STATUS in another column (filled with 0/1).
#' @return p-value of t-test for the two SURVIVAL_STATUS groups.
#' @export
assess_with_ttest <- function(X) {
    tt<-t.test(X[X$SURVIVAL_STATUS == 0, 1], X[X$SURVIVAL_STATUS == 1, 1])
    return(tt$p.value)
}


#' Make univariate GLM for assessing relevance of a variable 
#' 
#' @param X Data frame with numerical data in first column, and SURVIVAL_STATUS in another column (filled with 0/1).
#' @param plotting (optional). Whether to make a plot of the ROC curve. Default FALSE.
#' @return AUC of GLM for predicting the SURVIVAL_STATUS.
#' @export
assess_with_glm <- function(X, plotting=F) {
    model <- glm(SURVIVAL_STATUS ~., family=binomial(link='logit'), data=X)

    p <- predict(model, type="response")
    pr <- prediction(p, X$SURVIVAL_STATUS)

    prf <- performance(pr, measure = "tpr", x.measure = "fpr")
    if(plotting) {
        plot(prf)
    }

    auc <- performance(pr, measure = "auc")
    auc <- auc@y.values[[1]]
    return(auc)
}


#' Use ANOVA for assessing relevance of a variable to subtype
#' 
#' @param X Data frame with numerical data in first column, and CANCER_TYPE in another column.
#' @return p-value of ANOVA test.
#' @export
assess_subtype_with_anova <- function(X) {
    types_to_consider=c("ASTROCYTOMA","GBM","OLIGODENDROGLIOMA")
    Y <- X[X$CANCER_TYPE %in% types_to_consider,]
    colnames(Y) <- c("weights", "CANCER_TYPE")
    fitted <- aov( weights ~ CANCER_TYPE, data=Y)
    pvalue <- summary.aov(fitted)[[1]][["Pr(>F)"]][1]
    return(pvalue)
}


#' Use MGLM for assessing relevance of a variable to subtype
#' 
#' @param X Data frame with numerical data in first column, and CANCER_TYPE in another column.
#' @return accuracy value of fitted model.
#' @export
assess_subtype_with_mglm <- function(X) {
    types_to_consider=c("ASTROCYTOMA","GBM","OLIGODENDROGLIOMA")
    index = c(1:3)
    names(index) <- types_to_consider
    Y <- X[X$CANCER_TYPE %in% types_to_consider,]
    Y$CANCER_TYPE_numeric <- sapply(Y$CANCER_TYPE, FUN=function(x){index[x]})
    colnames(Y) <- c("weights", "CANCER_TYPE", "CANCER_TYPE_numeric")
    model <- multinom(CANCER_TYPE_numeric ~ weights, data=Y, trace = FALSE)
    pred <- predict(model)
    accuracy <- sum(pred==Y$CANCER_TYPE_numeric)/length(Y$CANCER_TYPE_numeric)
    return(accuracy)
}


#' Assess relevance of variables to outcome
#' 
#' @param df Data frame as returned by load_data. Rows are samples. Columns are features, including SURVIVAL_STATUS and CANCER_TYPE, as well as the input feature names.
#' @param feature_names Feature names (e.g. genes), as returned by load_data.
#' @param cores (optional) Number of cores to use for parallel computation. Default 1.
#' @return stats=list(pvalue_ttest, auc_glm, pvalue_anova, acc_mglm). These four elements are vectors with one entry for each feature (with names given by feature_names).
#' @export
calculate_stats <- function(df, feature_names, cores=1) {

    get_p_ttest <- function(name) { return(assess_with_ttest(         df[, c(name, "SURVIVAL_STATUS")] )) }
    get_auc <- function(name)     { return(assess_with_glm(           df[, c(name, "SURVIVAL_STATUS")] )) }
    get_p_anova <- function(name) { return(assess_subtype_with_anova( df[, c(name, "CANCER_TYPE")]  )) }
    get_acc <- function(name)     { return(assess_subtype_with_mglm(  df[, c(name, "CANCER_TYPE")]  )) }

    # pvalue_ttest <- unlist(pbmclapply(feature_names, FUN=get_p_ttest))
    # auc_glm      <- unlist(pbmclapply(feature_names, FUN=get_auc, mc.cores=cores))
    # pvalue_anova <- unlist(pbmclapply(feature_names, FUN=get_p_anova, mc.cores=cores))
    # acc_mglm     <- unlist(pbmclapply(feature_names, FUN=get_acc, mc.cores=cores))
    pvalue_ttest <- unlist(parallel::mclapply(feature_names, FUN=get_p_ttest))
    auc_glm      <- unlist(parallel::mclapply(feature_names, FUN=get_auc, mc.cores=cores))
    pvalue_anova <- unlist(parallel::mclapply(feature_names, FUN=get_p_anova, mc.cores=cores))
    acc_mglm     <- unlist(parallel::mclapply(feature_names, FUN=get_acc, mc.cores=cores))
    names(pvalue_ttest) <- feature_names
    names(auc_glm) <- feature_names
    names(pvalue_anova) <- feature_names
    names(acc_mglm) <- feature_names

    stats <- list(pvalue_ttest=pvalue_ttest, auc_glm=auc_glm, pvalue_anova=pvalue_anova, acc_mglm=acc_mglm)
    return(stats)
}


#' Sort features and write to file
#' 
#' Utility function
#' @param stats As returned by calculate_stats.
#' @param feature_names As returned by load_data.
#' @param write (optional) A string added to output filenames to aid in identification (e.g. write="ge", write="cn", or write="cnge"). Default NA, no output saved to file.
#' @return list(ttest_ordered, glm_ordered)
#' @export
rank_features <- function(stats, feature_names, write=NA) {
    ttest_ordered <- names(stats$pvalue_ttest)[order(stats$pvalue_ttest)]
    glm_ordered   <- names(stats$auc_glm)[order(stats$auc_glm, decreasing=T)]
    anova_ordered <- names(stats$pvalue_anova)[order(stats$pvalue_anova)]
    mglm_ordered  <- names(stats$acc_mglm)[order(stats$acc_mglm, decreasing=T)]

    if(!is.na(write)) {
        write.csv(ttest_ordered, paste0(write, "_ttest_ordered_features.csv"), quote=F, row.names=F)
        write.csv(glm_ordered, paste0(write, "_glm_ordered_features.csv"), quote=F, row.names=F)
        write.csv(anova_ordered, paste0(write, "_anova_ordered_features.csv"), quote=F, row.names=F)
        write.csv(mglm_ordered, paste0(write, "_mglm_ordered_features.csv"), quote=F, row.names=F)
    }
    ranked <- list(ttest_ordered=ttest_ordered,
                   glm_ordered=glm_ordered,
                   anova_ordered=anova_ordered,
                   mglm_ordered=mglm_ordered)
    return(ranked)
}


#' Combine rankings into one
#' 
#' More flexible than combining a fixed number from each.
#' @param rankings Data frame or matrix. Each column is a list of elements representing a ranking. The columns should have the same elements.
#' @param modality Either "intersection" or "union". Corresponds to the use of the max or min function on rankings for a given named element. The effect is the same as to form the intersection (respectively union) of the first K elements from each ranking, recording the order in which new elements are added.
#' @return A vector, the same length as the columns of ranking, representing the merged ranking.
#' @export
merge_rankings <- function(rankings, modality=c("intersection", "union")) {
    number_entries <- dim(rankings)[1]
    number_rankings <- dim(rankings)[2]
    feature_names <- sort(rankings[,1])
    for(i in 2:number_rankings) {
        if(!all(feature_names == sort(rankings[,i]))) {
            stop(paste0("Ranking 1 and ", i, " involve different named feature sets"))
        }
    }
    
    direct_ranks <- matrix(0, nrow=dim(rankings)[1], ncol=dim(rankings)[2])
    rownames(direct_ranks) <- feature_names
    for(j in 1:(dim(rankings)[2])) {
        temp_rank <- c(1:number_entries)
        names(temp_rank) <- rankings[,j]
        direct_ranks[,j] <- temp_rank[feature_names]
    }

    modality <- match.arg(modality)
    if(modality == "intersection") {
        summaryfunction <- max
    }
    if(modality == "union") {
        summaryfunction <- min
    }

    ranks_by_name <- sapply(1:number_entries, FUN=function(x){summaryfunction(direct_ranks[x,])})
    names(ranks_by_name) <- feature_names
    new_ranking <- names(sort(ranks_by_name))
    return(new_ranking)
}


#' Combine various selected rankings.
#' 
#' @param ttest_ordered As returned by stats.
#' @param glm_ordered As returned by stats.
#' @param anova_ordered As returned by stats.
#' @param mglm_ordered As returned by stats.
#' @param feature_names As returned by load_data.
#' @return Vector, merged ranking.
#' @export
merge_specific_rankings <- function(ttest_ordered, glm_ordered, anova_ordered, mglm_ordered, feature_names) {
    rankings_t_glm <- cbind(ttest_ordered, glm_ordered)
    merged_t_glm <- merge_rankings(rankings_t_glm, modality="intersection")

    rankings_anova_mglm <- cbind(anova_ordered, mglm_ordered)
    merged_anova_mglm <- merge_rankings(rankings_anova_mglm, modality="intersection")

    rankings_merged <- cbind(merged_t_glm, merged_anova_mglm)
    merged <- merge_rankings(rankings_merged, modality="union")
    
    return(merged)
}


#' Write results to file
#' 
#' Write subset of data involving selected features to file. Includes CANCER_TYPE and SURVIVAL_STATUS.
#' @param features Features to include.
#' @param df As returned by load_data.
#' @param write (optional) A string added to output filenames to aid in identification. Default NA, no output is saved to file.
#' @return The selected dataset.
#' @export
write_subset <- function(features, df, write=NA) {
    if(is.na(write)) {
        prepend <- ""
    } else {
        prepend <- paste0(write, "_")
    }
    filename <- paste0(prepend, "merged_ranked.csv")
    df_selected <- df[, c("SURVIVAL_STATUS",features)]
    if(!is.na(write)) {
        write.csv(df_selected, filename, quote=FALSE)
    }
    return(df_selected)
}


#' Make GLM model
#' 
#' @param X Data frame with with feature data and SURVIVAL_STATUS.
#' @param plotting (optional). Whether to make a plot of the ROC curve. Default FALSE.
#' @param quiet (optional). Suppress output to console. Default FALSE.
#' @return AUC of GLM for predicting the SURVIVAL_STATUS.
#' @export
build_glm <- function(X, plotting=FALSE, quiet=FALSE) {
    glm_model <- glm(SURVIVAL_STATUS ~., family=binomial(link='logit'), data=X)

    p <- predict(glm_model, type="response")
    pr <- prediction(p, X$SURVIVAL_STATUS)
    prf <- performance(pr, measure = "tpr", x.measure = "fpr")
    if(plotting) {
        plot(prf)
    }

    # fprs <- slot(prf, "x.values")[[1]]
    # tprs <- slot(prf, "y.values")[[1]]
    cutoffs <- slot(prf, "alpha.values")[[1]]
    accs <- sapply(cutoffs, FUN=function(cutoff){
                                sum((p > cutoff) == X$SURVIVAL_STATUS)/length(p)
                            })
    i <- which.max(accs)
    cutoff <- cutoffs[i]
    acc  <- sum((p > cutoff) == X$SURVIVAL_STATUS)/length(p)
    sens <- sum((p > cutoff) & X$SURVIVAL_STATUS)/sum(X$SURVIVAL_STATUS)
    spec <- sum(!(p > cutoff) & !X$SURVIVAL_STATUS)/sum(!X$SURVIVAL_STATUS)

    auc <- performance(pr, measure = "auc")
    auc <- auc@y.values[[1]]

    if(!quiet) {
        cat("Contingency table  :")
        print(table((p > cutoff) , X$SURVIVAL_STATUS ))
        cat("                *Actual survival status (0/1)\n")
        cat("                *Predicted death        (FALSE/TRUE)\n")
        cat("\n")

        cat(paste0("Accuracy           : ", round(acc, digits=4), "\n"))
        cat(paste0("Sensitivity        : ", round(sens, digits=4), "\n"))
        cat(paste0("Specificity        : ", round(spec, digits=4), "\n"))
        cat(paste0("AUC                : ", round(auc, digits=4), "\n"))
    }

    model_wrapper <- list()
    model_wrapper$cutoff <- cutoff
    model_wrapper$acc <- acc
    model_wrapper$sens <- sens
    model_wrapper$spec <- spec
    model_wrapper$auc <- auc
    model_wrapper$glm_model <- glm_model
    model_wrapper$coef <- summary(glm_model)$coef[,"Estimate"]
    return(model_wrapper)
}


#' Fold-creation for cross validation respecting subtype and survival status
#' 
#' @param df Data frame as returned by load_data or pipeline. Only needs columns CANCER_TYPE and SURVIVAL_STATUS to be present.
#' @param k (optional) Number of parts in the partition. Default 3. 
#' @return A partition (list) balanced with respect to CANCER_TYPE and SURVIVAL_STATUS
#' @export
create_folds_subtype <- function(df, k=3) {
    set.seed(123456)
    types_to_consider=c("ASTROCYTOMA","GBM","OLIGODENDROGLIOMA","UNKNOWN")
    fine_partitions <- lapply(types_to_consider, FUN=function(x){
                                    indices <- createFolds(df$SURVIVAL_STATUS[df$CANCER_TYPE == x], k=k)
                                    fine_partition <- lapply(indices, FUN=function(indices_part){ rownames(df)[df$CANCER_TYPE == x][indices_part] })
                                    return(fine_partition)
                                })
    partition <- mapply(c, fine_partitions[[1]], fine_partitions[[2]], fine_partitions[[3]], fine_partitions[[4]])
    return(partition)
}


#' Perform cross validation of whole pipeline
#' 
#' @param ld Data as returned by load_data or pipeline.
#' @param k (optional) Number of parts in the partition. Default 3. 
#' @param cores (optional) Number of cores for parallel computation. Default 1.
#' @param threshold (optional) Number of features to use. Default 100.
#' @param hierarchical_feature_selection (optional) Whether to use hierarchical clustering for final step of feature selection. Default FALSE.
#' @param number_hierarchical_features Number of final features to select using hierarchical clustering.
#' @return ...
#' @export
cross_validation <- function(ld, k=3, cores=1, threshold=100, hierarchical_feature_selection=FALSE, number_hierarchical_features) {
    df <- ld[[1]]
    feature_names <- ld[[2]]
    partition <- create_folds_subtype(df, k=k)
    cv <- lapply(partition, function(x){
                                training_fold <- df[setdiff(rownames(df),x), ]
                                test_fold <- df[x, ]
                                result <- pipeline( c(),
                                                    ld=list(training_fold, feature_names),
                                                    testset=test_fold,
                                                    drop.zeros=TRUE,
                                                    threshold=threshold,
                                                    cores=cores,
                                                    model_only=TRUE,
                                                    hierarchical_feature_selection=hierarchical_feature_selection,
                                                    number_hierarchical_features=number_hierarchical_features,
                                                    quiet=TRUE)
                                return(result)
                            })

    stats <- data.frame(matrix(0, nrow=4, ncol=k))
    for(j in 1:k) {
        stats[1,j] <- cv[[j]]$acc
        stats[2,j] <- cv[[j]]$sens
        stats[3,j] <- cv[[j]]$spec
        stats[4,j] <- cv[[j]]$auc
    }
    means <- apply(MARGIN=1, stats, mean)
    stats <- cbind(stats, means)
    colnames(stats) <- c(paste0("fold",c(1:k)), "mean")
    rownames(stats) <- c("acc","sens","spec","auc")
    cat("\n\n")
    cat("Cross-validation statistics (prediction on test data):\n")
    cat("\n")
    print(round(stats, digits=2))
    cat("\n")
    return(stats)
}


#' Select final representative features using hierarchical clustering
#' 
#' @param ld Data as returned by load_data or pipeline.
#' @param ranking Feature set ranking.
#' @param N Number of features to return.
#' @param M Number of features to start with.
#' @return Vector of N features.
#' @export
select_representative_features <- function(ld, ranking, N, M) {
    df <- ld[[1]]
    feature_names <- ld[[2]]
    df_selected <- df[, ranking[1:M]]

    temp_rank <- c(1:length(ranking))
    names(temp_rank) <- ranking
    direct_ranks <- temp_rank[feature_names]
    names(direct_ranks) <- feature_names

    clusters <- hclust(dist(t(df_selected), method = "euclidean"), method = "complete", members = NULL)
    cluster_cut <- cutree(clusters, N)
    representatives <- sapply(unique(cluster_cut), FUN=function(j){
                                    members <- names(which(cluster_cut == j))
                                    representative <- names(which.min(direct_ranks[members]))
                                    return(representative)
                                })
    return(representatives[order(direct_ranks[representatives])])
}


#' Pipeline
#' 
#' @param filenames Filenames of the phenotype, outcome, and feature matrix (for challenges 1/2/3).
#' @param ld (optional). Supply the loaded data (as returned by a prior pipeline(...) call), if desired.
#' @param merged_ranking (optional). Supply the computed merged_ranking, to skip the pipeline forward to the last steps.
#' @param plotting (optional). Whether to make a plot of the ROC curves. Default FALSE.
#' @param drop.zeros Whether to drop features that are identically zero.
#' @param write (optional) A string added to output filenames to aid in identification. Default NA, no addition.
#' @param threshold Number of features to use in ranking stage.
#' @param cores Number of cores for parallel computation.
#' @param model_only (optional) Whether to return only the model/stats. Default FALSE.
#' @param hierarchical_feature_selection (optional) Whether to use hierarchical clustering for final step of feature selection. Default FALSE.
#' @param number_hierarchical_features Number of final features to select using hierarchical clustering.
#' @param testset (optional). If provided, the reporting on the final model is for its performance on the testset. Default NULL.
#' @param quiet (optional). Suppress outpute to console. Default FALSE.
#' @return ld (see load_data). If model_only=TRUE, returns only summary statistics for the run (for validation, etc.).
#' @export
pipeline <- function(filenames, ld=NULL, merged_ranking=NULL, plotting=FALSE, drop.zeros=FALSE, write=NA, threshold=100, cores=1, model_only=FALSE, hierarchical_feature_selection=FALSE, number_hierarchical_features=NA, testset=NULL, quiet=FALSE) {
    #Load data
    if(is.null(ld)) {
        ld <- load_data(filenames, drop.zeros=drop.zeros)
    }
    df <- ld[[1]]
    feature_names <- ld[[2]]

    #Feature assessment
    if(is.null(merged_ranking)) {
        stats <- calculate_stats(df, feature_names, cores=cores) 
        ranked <- rank_features(stats, feature_names, write=write)
        merged_ranking <- merge_specific_rankings(ranked$ttest_ordered, ranked$glm_ordered, ranked$anova_ordered, ranked$mglm_ordered,feature_names)
        ld[[3]] <- merged_ranking
        if(!is.na(write)) {
            write.table(merged_ranking, paste0(write, "_all_features_ranked.csv"), col.names=FALSE, row.names=FALSE, quote=FALSE)
        }
    }

    #Feature selection
    if(hierarchical_feature_selection) {
        if(is.na(number_hierarchical_features)) {
            stop("If using hierarchical_feature_selection, need to supply number_hierarchical_features.")
        }
        number_of_features <- number_hierarchical_features
        final_features <- select_representative_features(ld, merged_ranking, number_of_features, threshold)
    } else {
        number_of_features <- threshold
        final_features <- merged_ranking[1:number_of_features]
    }
    
    #Build model
    if(!is.na(write)) {
        write.table(final_features, paste0(write, "_features.csv"), col.names=FALSE, row.names=FALSE, quote=FALSE)
        write_subset(final_features, df, write=write)
    }
    if(!quiet) {
        cat("\n")
        cat(paste0("GLM model for ", filenames[1], ", ... (", write, ")", "\n\n"))
        
        cat(paste0("Number of features : ", number_of_features))
        if(hierarchical_feature_selection) {
            cat(paste0("        (Hierarchical cluster representatives)"))
        }
        cat("\n")
        cat(paste(final_features[1:(min(length(final_features), 10))]))
        cat("... \n")
    }

    model_wrapper <- build_glm(df[, c("SURVIVAL_STATUS", "CANCER_TYPE", final_features)], plotting=plotting, quiet=quiet)
    needed_data <- c(model_wrapper$cutoff, model_wrapper$coef)
    names(needed_data)[1] <- "(Cutoff)"
    if(!is.na(write)) {
        file <- paste0(write, "_coefficients")
        write.table(needed_data, paste0(file, ".tsv"), col.names=FALSE, sep="\t", quote=FALSE)
        if(!quiet) {
            cat("\n")
            cat(paste0("Wrote data needed to run model to ", paste0(file, ".tsv"), "\n"))
        }
        # zip(paste0(file, ".zip"), paste0(file, ".tsv"))
        # cat(paste0("Compressed file is ", paste0(file, ".zip"), "\n"))
    }

    if(!is.null(testset)) {
        p <- predict(model_wrapper$glm_model, newdata=testset, type="response")
        pr <- prediction(p, testset$SURVIVAL_STATUS)

        prf <- performance(pr, measure = "tpr", x.measure = "fpr")
        if(plotting) {
            plot(prf)
        }

        # fprs <- slot(prf, "x.values")[[1]]
        # tprs <- slot(prf, "y.values")[[1]]
        cutoffs <- slot(prf, "alpha.values")[[1]]
        accs <- sapply(cutoffs, FUN=function(cutoff){
                                    sum((p > cutoff) == testset$SURVIVAL_STATUS)/length(p)
                                })
        i<-which.max(accs)
        cutoff <- cutoffs[i]
        acc  <- sum((p > cutoff) == testset$SURVIVAL_STATUS)/length(p)
        sens <- sum((p > cutoff) & testset$SURVIVAL_STATUS)/sum(testset$SURVIVAL_STATUS)
        spec <- sum(!(p > cutoff) & !testset$SURVIVAL_STATUS)/sum(!testset$SURVIVAL_STATUS)

        auc <- performance(pr, measure = "auc")
        auc <- auc@y.values[[1]]

        predict_stats <- list()
        predict_stats$acc <- acc
        predict_stats$sens <- sens
        predict_stats$spec <- spec
        predict_stats$auc <- auc
        return(predict_stats)
    }

    if(model_only) {
        return(model_wrapper)
    } else {
        return(ld)
    }
}







