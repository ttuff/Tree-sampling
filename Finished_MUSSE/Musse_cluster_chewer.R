library(ape)
library(apTreeshape )
library(phytools)
library(diversitree)
library(Rcpp)
library(picante)
library(diversitree)

## sample file ###############################
pruned_trees <- list.files(path="/scratch/ttuff/Tree_sampling/Finished_tree_prunes")
finished_musse <- list.files(path="/scratch/ttuff/Tree_sampling/Finished_MUSSE_sets")
already_done <- which(!(pruned_trees %in% finished_musse))

if(length(already_done) != 0)(unique_remaining <- pruned_trees[already_done])
pre_files <-  sample(unique_remaining, 1)

files <- paste0("/scratch/ttuff/Tree_sampling/Finished_tree_prunes/", pre_files)
load(as.character(files))

stringer <- strsplit(pre_files, "_")[[1]]
phy <- pruned_trees_and_stats[[2]]
FRACTION_REMAINING <- 48
## Run Musse ###############################


# phy <- working_tree 
chew_musse <- function(phy, FRACTION_REMAINING = 1 , NUM_STATES=4, STEP_NUM=1000){
  print("just entered chew_musse")
  flag<-runif(1)
	print("just ran runif(1)")
	# starting point parameter estimates for MLE and MCMC.
	startingPoint <- starting.point.musse(phy, NUM_STATES)
	# create a prior distribution for MCMC
	exp_prior <- make.prior.exponential(1/ (2*startingPoint[1]))
	print("made starting point")
	sprintf("%s : made startingPoint and prior", flag)
	
	# DO analysis with sampling correction.
	# make likelihood function accounting for fraction sampled
	lik.pruneAdjusted <- make.musse (phy, phy$tip.state, NUM_STATES,
		sampling.f=rep(FRACTION_REMAINING/100, NUM_STATES)
	)
	# manually constrain the model to realize its multiple binary traits.
	lik.pruneAdjusted.multitrait <- constrain(lik.pruneAdjusted,  q14 ~ 0, q23 ~ 0, q32 ~ 0, q41 ~ 0)
	# Find MLE answer
	MLE_answer.pruneAdjusted<- find.mle(lik.pruneAdjusted.multitrait,
		startingPoint[argnames(lik.pruneAdjusted.multitrait)])
	# do MCMC chain.
	sprintf("%s : found MLE for adjusted analysis", flag)
	samples.pruneAdjusted <- mcmc(lik.pruneAdjusted.multitrait, 
		coef(MLE_answer.pruneAdjusted), 
		nsteps=STEP_NUM, 
		prior = exp_prior, 
		w=.1, 
		print.every=0
	)
	sprintf("%s : finished adjusted mcmc", flag)
	
	# DO analysis WITHOUT correction.
	# make likelihood function
	lik.notAd <- make.musse (phy, phy$tip.state, NUM_STATES)
	# manually constrain the model to realize its multiple binary traits.
	lik.notAd.multitrait <- constrain(lik.notAd,  q14 ~ 0, q23 ~ 0, q32 ~ 0, q41 ~ 0)
	# Find MLE answer
	sprintf("%s : starting non-adjusted MLE", flag)
	MLE_answer.notAd<- find.mle(lik.notAd.multitrait,
		startingPoint[argnames(lik.notAd.multitrait)])
	# do MCMC chain.
	sprintf("%s : finished non-adjusted MLE", flag)
	samples.notAd <- mcmc(lik.notAd.multitrait, 
		coef(MLE_answer.notAd), 
		nsteps=STEP_NUM, 
		prior = exp_prior, 
		w=.1, 
		print.every=0
	)
	sprintf("%s : finished non-adjusted mcmc, started stats", flag)
	# handles errors.
	
	sprintf("%s : finished stats", flag)
	
	return(list( MLE_answer.notAd, samples.notAd, MLE_answer.pruneAdjusted, samples.pruneAdjusted))
}





print("started")
fraction_removed <- as.numeric(stringer[2])
try(musse_analysis_1a <- chew_musse(pruned_trees_and_stats[[2]], FRACTION_REMAINING = fraction_removed), silent=TRUE)
print("1_done")
try(musse_analysis_2a <- chew_musse(pruned_trees_and_stats[[3]], FRACTION_REMAINING = fraction_removed), silent=TRUE)
print("2_done")
try(musse_analysis_3a <- chew_musse(pruned_trees_and_stats[[4]], FRACTION_REMAINING = fraction_removed), silent=TRUE)
print("3_done")
try(musse_analysis_4a <- chew_musse(pruned_trees_and_stats[[5]], FRACTION_REMAINING = fraction_removed), silent=TRUE)
print("4_done")
try(musse_analysis_5a <- chew_musse(pruned_trees_and_stats[[6]], FRACTION_REMAINING = fraction_removed), silent=TRUE)
print("5_done")
try(musse_analysis_6a <- chew_musse(pruned_trees_and_stats[[7]], FRACTION_REMAINING = fraction_removed), silent=TRUE)
print("6_done")

musse_analysis_1 <- musse_analysis_1a[[1]]$par 
musse_analysis_1_corrected <- musse_analysis_1a[[3]]$par
musse_analysis_2 <- musse_analysis_2a[[1]]$par 
musse_analysis_2_corrected <- musse_analysis_2a[[3]]$par
musse_analysis_3 <- musse_analysis_3a[[1]]$par
musse_analysis_3_corrected <- musse_analysis_3a[[3]]$par
musse_analysis_4 <- musse_analysis_4a[[1]]$par 
musse_analysis_4_corrected <- musse_analysis_4a[[3]]$par
musse_analysis_5 <- musse_analysis_5a[[1]]$par 
musse_analysis_5_corrected <- musse_analysis_5a[[3]]$par
musse_analysis_6 <- musse_analysis_6a[[1]]$par
musse_analysis_6_corrected <- musse_analysis_6a[[3]]$par

musse_analysis_matrix <- rbind(musse_analysis_1, musse_analysis_2, musse_analysis_3, musse_analysis_4, musse_analysis_5, musse_analysis_6)

musse_analysis_matrix_corrected <- rbind(musse_analysis_1_corrected, musse_analysis_2_corrected, musse_analysis_3_corrected, musse_analysis_4_corrected, musse_analysis_5_corrected, musse_analysis_6_corrected)


## build summary matrix ###############################
summary_matrix <- matrix(rep(NA, 75*6), 6, 75)
colnames(summary_matrix)  <- c("Prune_type", "Percent_removed", "Original_tree_size", "lambda1_input", "lambda2_input", "lambda3_input", "lambda4_input", "mu1_input", "mu2_input", "mu3_input", "mu4_input", "q12_input", "q13_input", "q14_input","q21_input",  "q23_input","q24_input", "q31_input","q32_input",  "q34_input",  "q41_input", "q42_input","q43_input", "lambda1_output_uncorrected", "lambda2_output_uncorrected", "lambda3_output_uncorrected", "lambda4_output_uncorrected", "mu1_output_uncorrected", "mu2_output_uncorrected", "mu3_output_uncorrected", "mu4_output_uncorrected", "q12_output_uncorrected", "q13_output_uncorrected",  "q21_output_uncorrected",  "q24_output_uncorrected", "q31_output_uncorrected",  "q34_output_uncorrected",  "q42_output_uncorrected", "q43_output_uncorrected",  "lambda1_output_corrected", "lambda2_output_corrected", "lambda3_output_corrected", "lambda4_output_corrected", "mu1_output_corrected", "mu2_output_corrected", "mu3_output_corrected", "mu4_output_corrected", "q12_output_corrected", "q13_output_corrected",  "q21_output_corrected",  "q24_output_corrected", "q31_output_corrected",  "q34_output_corrected",  "q42_output_corrected", "q43_output_corrected",  "number_of_branches", "Pylo_diversity_is_sum_of_BL", "average_phylogenetic_diversity_is_mean_of_BL", "variance_Pylo_diversity_is_variance_of_BL", "F_quadratic_entropy_is_sum_of_PD", "Mean_pairwise_distance", "variance_pairwise_distance", "colless_stat", "sackin_index", "tree_shape_stat", "Evolutionary_distinctiveness_sum", "mean_Phylogenetic_isolation", "variance_Phylogenetic_isolation","gamma", "gamma_p_value", "speciation_rate", "extinction_rate", "extinction_per_speciation", "speciation_minus_extinction", "filename"
)


statistics <- as.matrix(pruned_trees_and_stats[[8]])
dim(statistics)
summary_matrix[,1] <- c("full_tree", "random_prune", "shallow_trimmed_prune","deep_trimmed_prune", "overdispersed_prune", "underdispersed_prune")
summary_matrix[,2] <- stringer[2]
summary_matrix[,3] <- stringer[6]

summary_matrix[,4] <- stringer[8]
summary_matrix[,5] <- stringer[10]
summary_matrix[,6] <- stringer[12]
summary_matrix[,7] <- stringer[14]
summary_matrix[,8] <- stringer[16]
summary_matrix[,9] <- stringer[18]
summary_matrix[,10] <- stringer[20]
summary_matrix[,11] <- stringer[22]
summary_matrix[,12] <- stringer[24]
summary_matrix[,13] <- stringer[26]
summary_matrix[,14] <- stringer[28]
summary_matrix[,15] <- stringer[30]
summary_matrix[,16] <- stringer[32]
summary_matrix[,17] <- stringer[34]
summary_matrix[,18] <- stringer[36]
summary_matrix[,19] <- stringer[38]
summary_matrix[,20] <- stringer[40]
summary_matrix[,21] <- stringer[42]
summary_matrix[,22] <- stringer[44]
summary_matrix[,23] <- stringer[46]
#summary_matrix[,24] <- stringer[48]
#load(as.character(available_files[i]))
summary_matrix[1:6,24:39] <- musse_analysis_matrix
summary_matrix[1:6,40:55] <- musse_analysis_matrix_corrected
summary_matrix[1:6,56:74] <- statistics[,2:20]
summary_matrix[,75] <- as.character(pre_files)

## save ###############################
pruned_trees_and_stats_and_musse <- list(pruned_trees_and_stats[[2]], pruned_trees_and_stats[[3]], pruned_trees_and_stats[[4]], pruned_trees_and_stats[[5]], pruned_trees_and_stats[[6]], pruned_trees_and_stats[[7]], summary_matrix, pruned_trees_and_stats[[9]], pruned_trees_and_stats[[10]])

save(pruned_trees_and_stats_and_musse, file=paste("/scratch/ttuff/Tree_sampling/Finished_MUSSE_sets/", pre_files, sep=""))

