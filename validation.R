library(bbutilities)  # Uses load_data
options(digits=4)

simplify_type <- function(cancer_type) {
    if(cancer_type %in% c("ASTROCYTOMA", "GBM", "OLIGODENDROGLIOMA", "UNKNOWN_MIXED_UNCLASSIFIED")) {
        return(cancer_type)
    }
    warning(paste0("Unexpected CANCER_TYPE value: ", cancer_type))
    return("UNKNOWN_MIXED_UNCLASSIFIED")
}

evaluate_model_from_coefficients <- function(coefficients, data, challenge=0) {
    cutoffs_row <- which(rownames(coefficients) == "(Cutoff)")
    intercepts_row <- which(rownames(coefficients) == "(Intercept)")

    cutoffs <- coefficients[cutoffs_row,]
    intercepts <- coefficients[intercepts_row,]

    model_mat <- coefficients[-c(cutoffs_row, intercepts_row),]
    features <- rownames(model_mat)
    if(challenge == 2) {
        features <- c(paste0("X", features))
    }
    rownames(model_mat) <- features
    feature_mat <- data[, c(features, "CANCER_TYPE"),]
    predictions <- apply(feature_mat, MARGIN=1, function(x) {
            cancer_type <- simplify_type(x["CANCER_TYPE"])
            cutoff <- cutoffs[cancer_type]
            intercept <- intercepts[cancer_type]
            value <- intercept + sum(sapply(features, function(f) { as.numeric(x[f])*model_mat[f, cancer_type] } ))
            return(value > cutoff)
        })
    return(predictions)
}

evaluate_from_R_glm <- function(glm_models, cutoffs, data) {
    cancer_types <- data$CANCER_TYPE
    if(!all(sort(unique(names(glm_models))) == sort(unique(cancer_types)))) {
        warning("New data CANCER_TYPE possible values don't match expected values.")
    }
    
    predictions_by_type <- lapply(names(glm_models), function(cancer_type) {
            glm_model <- glm_models[[cancer_type]]
            cutoff <- cutoffs[[cancer_type]]
            newdata <- data[data$CANCER_TYPE == cancer_type, ]
            p <- predict(glm_model, newdata=newdata, type="response")
            return(sapply(p > cutoff, function(bool){if(bool) {1} else {0}} ))
        })
    predictions <- do.call(c, predictions_by_type)
    if(!all(sort(unique(rownames(data))) == sort(unique(names(predictions))))) {
        warning("Resulting predictions have wrong patient IDs.")
    }
    predictions <- predictions[rownames(data)]
    return(predictions)
}

rds_files <- paste0(paste0("subchallenge", c(1,2,3)), "_glm_models.rds")
rds_data <- lapply(rds_files, function(file){ readRDS(file) })
glm_models <- lapply(c(1:3), function(j) { rds_data[[j]][[1]] })
cutoffs <- lapply(c(1:3), function(j) { rds_data[[j]][[2]] })

check_models_on_training_data <- function(glm_models, cutoffs) {
    f1 <- c("sc1_Phase1_GE_Phenotype.tsv", "sc1_Phase1_GE_Outcome.tsv", "sc1_Phase1_GE_FeatureMatrix.tsv")
    f2 <- c("sc2_Phase1_CN_Phenotype.tsv", "sc2_Phase1_CN_Outcome.tsv", "sc2_Phase1_CN_FeatureMatrix.tsv")
    f3 <- c("sc3_Phase1_CN_GE_Phenotype.tsv", "sc3_Phase1_CN_GE_Outcome.tsv", "sc3_Phase1_CN_GE_FeatureMatrix.tsv")
    filename_groups <- list(f1, f2, f3)
    data_frames <- lapply(filename_groups, function(filenames) {
            ld <- load_data(filenames)
            return(ld[[1]])
        })

    coefficient_files <- paste0(paste0("subchallenge", c(1,2,3)), "_coefficients.tsv")
    model_coefficients <- lapply(coefficient_files, function(file){ read.table(file, sep="\t") })

    predictions <- lapply(c(1:3), function(j) {
            # evaluate_model_from_coefficients(model_coefficients[[j]], data_frames[[j]], challenge=j)
            evaluate_from_R_glm(glm_models[[j]], cutoffs[[j]], data_frames[[j]])
        })

    labels <- lapply(data_frames, function(df) {
            survival <- df[, "SURVIVAL_STATUS"]
            names(survival) <- rownames(df)
            return(survival)
        })

    combined <- lapply(c(1:3), function(j) {
            cbind(predictions[[j]], labels[[j]])
        })

    for(j in 1:3) {
        print(table(combined[[j]][,1], combined[[j]][,2]))
    }
    compare_predictions_ground_truth <- combined
    return(compare_predictions_ground_truth)
}

evaluate_models_on_test_data <- function(glm_models, cutoffs) {
    g1 <- c("sc1_Phase2_GE_Phenotype.tsv", "", "sc1_Phase2_GE_FeatureMatrix.tsv")
    g2 <- c("sc2_Phase2_CN_Phenotype.tsv", "", "sc2_Phase2_CN_FeatureMatrix.tsv")
    g3 <- c("sc3_Phase2_CN_GE_Phenotype.tsv", "", "sc3_Phase2_CN_GE_FeatureMatrix.tsv")
    filename_groups_phase2 <- list(g1, g2, g3)
    data_frames_phase2 <- lapply(filename_groups_phase2, function(filenames) {
            ld <- load_data(filenames, has_labels=FALSE)
            return(ld[[1]])
        })

    predictions_phase2 <- lapply(c(1:3), function(j) {
            evaluate_from_R_glm(glm_models[[j]], cutoffs[[j]], data_frames_phase2[[j]])
        })

    for(j in c(1:3)) {
        print(table(predictions_phase2[[j]]))
    }

    prediction_filenames <- paste0("subchallenge_", c(1,2,3), ".tsv")
    for(j in c(1:3)) {
        p <- predictions_phase2[[j]]
        write.table(cbind(PATIENT_ID=names(p), SURVIVAL_STATUS=p), prediction_filenames[j], sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    }
    return(predictions_phase2)
}


compare_predictions_ground_truth <- check_models_on_training_data(glm_models, cutoffs)
predictions_phase2 <- evaluate_models_on_test_data(glm_models, cutoffs)


