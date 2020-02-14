

D <- data.frame(matrix(0, nrow=0, ncol=7))

for(ch in 1:3) {
	for(st in c("ASTROCYTOMA", "GBM", "OLIGODENDROGLIOMA", "UNKNOWN_MIXED_UNCLASSIFIED")) {
		names1 <- c("(Cutoff)",names(glm_models[[ch]][[st]]$coefficients))
		cs1 <- c(cutoffs[[ch]][[st]], glm_models[[ch]][[st]]$coefficients)
		names2 <- rownames(model_coefficients[[ch]])
		cs2 <- model_coefficients[[ch]][,st]
		tt<-cbind(
			names1,
			names2,
			names1 == names2,
			cs1,
			cs2,
			cs1-cs2,
			abs(cs1-cs2) < 0.0000000001
		)
		D <- rbind(D, tt)
	}
}

colnames(D) <- c()

# The coefficients match. The feature names match (modulo X).
