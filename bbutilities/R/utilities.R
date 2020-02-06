
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
#' @import crayon

library(stats)
library(utils)
library(pbmcapply)
library(parallel)
library(ROCR)
library(nnet)
library(caret)
library(crayon)

#' Load data (from brain biomarker challenge files)
#'
#' @param filenames Vector of filenames: phenotype, outcome, and feature data.
#' @param drop.zeros Whether to drop features that are identically zero.
#' @param merge_extraneous Whether to merge UNKNOWN, MIXED, and UNCLASSIFIED into one category.
#' @return list(df, feature_names). df$type, df$outcome, df$gene1, df$gene2, etc.
#' @export
load_data <- function(filenames, drop.zeros=FALSE, merge_extraneous=TRUE) {
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

    if(merge_extraneous) {
        extra <- c("UNKNOWN", "MIXED", "UNCLASSIFIED")
        newlabel <- paste(extra, collapse="_")
        df[df$CANCER_TYPE %in% extra,"CANCER_TYPE"] <- newlabel
    }

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

    rankings_merged <- cbind(merged_anova_mglm, merged_t_glm)
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
    suppressWarnings(
    glm_model <- glm(SURVIVAL_STATUS ~., family=binomial(link='logit'), data=X)
    )

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
        # cat("Contingency table  :")
        # print(table((p > cutoff) , X$SURVIVAL_STATUS ))
        # # cat(paste0(italic("                *Actual survival status (0/1)"),"\n"))
        # # cat(paste0(italic("                *Predicted death        (FALSE/TRUE)"),"\n"))
        # cat("\n")

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
    model_wrapper$predictions <- data.frame(id=rownames(X), predictions=(p > cutoff), SURVIVAL_STATUS=X$SURVIVAL_STATUS)
    return(model_wrapper)
}


#' Perform cross validation of simple GLM for SURVIVAL_STATUS, from a given feature set 
#' 
#' @param features The features to use. 
#' @param df The data frame.
#' @param cancer_type (optional) Cancer type or types to restrict the model to. Default NULL.
#' @param k (optional) Number of folds for cross validation. Default 3.
#' @param seed (optional) Seed number for pseudo-random reproducibility. Default 12345.
#' @return A data frame with the accuracy, sensitivity, specificity, and AUC of the test predictions. Also includes mean summary over the folds.
#' @export
cross_validation_glm <- function(features, df, cancer_type=NULL, k=3, seed=12345) {
    df_small <- df[,c("SURVIVAL_STATUS",features)]
    if(!is.null(cancer_type)) {
        df_small <- df_small[df$CANCER_TYPE %in% cancer_type,]
    }
    set.seed(seed)
    partition <- createFolds(df_small$SURVIVAL_STATUS, k=k)
    N <- length(rownames(df_small))

    test_part <- function(part) {
        training_fold <- df_small[setdiff(1:N, part), ]
        test_fold <- df_small[part, ]
        suppressWarnings(
        model_wrapper <- build_glm(training_fold, plotting=FALSE, quiet=TRUE)
        )
        suppressWarnings(
        p <- predict(model_wrapper$glm_model, newdata=test_fold, type="response")
        )
        pr <- prediction(p, test_fold$SURVIVAL_STATUS)
        prf <- performance(pr, measure = "tpr", x.measure = "fpr")
        auc <- performance(pr, measure = "auc")
        auc <- auc@y.values[[1]]

        X <- test_fold
        cutoffs <- slot(prf, "alpha.values")[[1]]
        accs <- sapply(cutoffs, FUN=function(cutoff){
                                    sum((p > cutoff) == X$SURVIVAL_STATUS)/length(p)
                                })
        i <- which.max(accs)
        cutoff <- cutoffs[i]
        acc  <- sum((p > cutoff) == X$SURVIVAL_STATUS)/length(p)
        sens <- sum((p > cutoff) & X$SURVIVAL_STATUS)/sum(X$SURVIVAL_STATUS)
        spec <- sum(!(p > cutoff) & !X$SURVIVAL_STATUS)/sum(!X$SURVIVAL_STATUS)

        return(c(Accuracy=acc, Sensitivity=sens, Specificity=spec, AUC=auc))
    }

    cv_result <- data.frame(lapply(partition, test_part))
    cv_result <- cbind(cv_result, sapply(1:4, function(j){round( mean(as.numeric(cv_result[j,])), digits=3)}) )
    colnames(cv_result)[4] <- "means"
    return(cv_result)
}


#' Save model parameters
#' 
#' @param model_wrapper A GLM as returned by build_glm.
#' @param write An identifying string for output file name.
#' @export
write_model_parameters_to_file <- function(model_wrapper, write) {
    needed_data <- c(model_wrapper$cutoff, model_wrapper$coef)
    names(needed_data)[1] <- "(Cutoff)"
    if(!is.na(write)) {
        file <- paste0(write, "_coefficients")
        write.table(needed_data, paste0(file, ".tsv"), col.names=FALSE, sep="\t", quote=FALSE)
    }
}


sanitize_R_prepends <- function(features) {
    if(length(grep("^X", features)) == length(features)) {
        return(sub("^X", "", features))
    }
    else {
        return(features)
    }
}


#' Build final model and report statistics
#' 
#' @param df The data frame.
#' @param final_features Features to use for final model.
#' @param write An identifying string for output file names.
#' @param filenames For reporting.
#' @return cv_mean_auc Mean AUC among cross-validations in each submodel.
#' @export
build_and_report_final_model <- function(df, final_features, write, filenames) {
    cat("\n")
    cat(paste0(bold("GLM model for ", filenames[1], ", ... (", write, ")"), "\n\n"))

    types <- unique(df$CANCER_TYPE)
    types <- setdiff(types, c("UNCLASSIFIED", "MIXED"))

    model_wrappers <- lapply(types, function(type) {
            cat("\n")
            cat(paste0(bold$bgGreen(type), "\n"))
            model_wrapper <- build_glm(df[df$CANCER_TYPE == type, c("SURVIVAL_STATUS", final_features)], plotting=FALSE)
            cat("\n")
            cat("Cross validation\n")
            cv_results <- cross_validation_glm(final_features, df, cancer_type=type)
            print(cv_results)
            write_model_parameters_to_file(model_wrapper, paste0(write, "_", type))
            model_wrapper$cv_mean_auc <- cv_results$means[4]
            return(model_wrapper)
        })

    all_predictions <- data.frame(matrix(0, nrow=0, ncol=3))
    for(type in types) {
        all_predictions <- rbind(all_predictions, model_wrappers[[which(types == type)]]$predictions)
    }

    pred <- all_predictions$predictions
    survival <- all_predictions$SURVIVAL_STATUS
    acc  <- sum(pred == survival)/length(pred)
    sens <- sum(pred & survival)/sum(survival)
    spec <- sum(!pred & !survival)/sum(!survival)

    cat("\n")
    cat(paste0(bold$bgCyan("Full model stats"),"\n"))
    cat("Contingency table  :")
    print(table((pred) , survival))
    cat(paste0(italic("                *Actual survival status (0/1)"),"\n"))
    cat(paste0(italic("                *Predicted death        (FALSE/TRUE)"),"\n"))
    cat("\n")

    cat(paste0("Accuracy              : ", round(acc, digits=4), "\n"))
    cat(paste0("Sensitivity           : ", round(sens, digits=4), "\n"))
    cat(paste0("Specificity           : ", round(spec, digits=4), "\n"))
    mean_auc <- mean(sapply(1:length(model_wrappers), function(j){return(model_wrappers[[j]]$auc)} ))
    cat(paste0("Mean AUC of submodels : ", round(mean_auc, digits=4), "\n"))
    cv_mean_auc <- mean(sapply(1:length(model_wrappers), function(j){return(model_wrappers[[j]]$cv_mean_auc)} ))
    cat(paste0("Mean AUC (submodel CV): ", round(cv_mean_auc, digits=4), "\n"))

    cat("\n")
    cat(paste0(bold("Features:"),"\n"))
    cat(paste(sanitize_R_prepends(final_features), collapse="\n"))
    cat("\n\n")
    return(c(cv_mean_auc=cv_mean_auc, acc=acc))
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
#' @return ld (see load_data).
#' @export
pipeline <- function(filenames, ld=NULL, merged_ranking=NULL, plotting=FALSE, drop.zeros=FALSE, write=NA, threshold=100, cores=1) {
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
    number_of_features <- threshold
    final_features <- merged_ranking[1:number_of_features]

    #Build model
    acc_and_cv_mean_auc <- build_and_report_final_model(df, final_features, write, filenames)
    ld[[4]] <- acc_and_cv_mean_auc
    return(ld)
}




