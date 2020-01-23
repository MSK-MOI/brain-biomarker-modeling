library(bbutilities)

# Challenge 1
filenames <- c("sc1_Phase1_GE_Phenotype.tsv", "sc1_Phase1_GE_Outcome.tsv", "sc1_Phase1_GE_FeatureMatrix.tsv")
ld <- pipeline(filenames, threshold=50, drop.zeros=TRUE, cores=3, write="subchallenge1", hierarchical_feature_selection=FALSE)

# Challenge 2
filenames <- c("sc2_Phase1_CN_Phenotype.tsv", "sc2_Phase1_CN_Outcome.tsv", "sc2_Phase1_CN_FeatureMatrix.tsv")
ld <- pipeline(filenames, threshold=50, drop.zeros=TRUE, cores=3, write="subchallenge2", hierarchical_feature_selection=FALSE)

# Challenge 3
filenames <- c("sc3_Phase1_CN_GE_Phenotype.tsv", "sc3_Phase1_CN_GE_Outcome.tsv", "sc3_Phase1_CN_GE_FeatureMatrix.tsv")
ld <- pipeline(filenames, threshold=20, drop.zeros=TRUE, cores=3, write="subchallenge3", hierarchical_feature_selection=FALSE)


