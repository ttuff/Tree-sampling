library(diversitree)
library(ape)
library(apTreeshape )
library(phytools)
library(diversitree)
library(Rcpp)
library(picante)



##Functions




###### Specify function ##############################

# BirthDeath function
bd <- function(tree) {
  tree$edge.length <- tree$edge.length / max(tree$edge.length) 
  x <- birthdeath(tree)  
  b <- x$para[2] / (1 - x$para[1])
  d <- b - x$para[2]
  c(setNames(c(b, d), c("b", "d")), x$para)
}





evol.distinct2 <- function(tree, type = c("equal.splits", "fair.proportion"), 
                           scale = FALSE, use.branch.lengths = TRUE) 
{
  type <- match.arg(type)
  if (is.rooted(tree) == FALSE) 
    warning("A rooted phylogeny is required for meaningful output of this function", 
            call. = FALSE)
  if (scale == TRUE) {
    if (is.ultrametric(tree) == TRUE) 
      tree$edge.length <- tree$edge.length/(as.numeric(branching.times(tree)[1]))
    else tree$edge.length <- tree$edge.length/sum(tree$edge.length)
  }
  if (use.branch.lengths == FALSE) 
    tree$edge.length <- rep(1, length(tree$edge.length))
  for (i in 1:length(tree$tip.label)) {
    spp <- tree$tip.label[i]
    nodes <- .get.nodes(tree, spp)
    nodes <- nodes[1:(length(nodes) - 1)]
    internal.brlen <- tree$edge.length[which(tree$edge[, 
                                                       2] %in% nodes)]
    if (length(internal.brlen) != 0) {
      internal.brlen <- internal.brlen * switch(type, equal.splits = sort(rep(0.5, 
                                                                              length(internal.brlen))^c(1:length(internal.brlen))), 
                                                fair.proportion = {
                                                  for (j in 1:length(nodes)) {
                                                    sons <- .node.desc(tree, nodes[j])
                                                    n.descendents <- length(sons$tips)
                                                    if (j == 1) portion <- n.descendents else portion <- c(n.descendents, 
                                                                                                           portion)
                                                  }
                                                  1/portion
                                                })
    }
    ED <- sum(internal.brlen, tree$edge.length[which.edge(tree, 
                                                          spp)])
    if (i == 1) 
      w <- ED
    else w <- c(w, ED)
  }
  return(w)
}



# Phylogenetic signal
Dsig <- function(mytree, myWorld) {
  traits <- data.frame("trait" = myWorld[, 6],
                       "tips" = myWorld[, 8])
  compdata <- comparative.data(mytree, traits, 'tips')
  # Phylogenetic signal for binary traits (D of Fritz and Purvis 2010)
  trait <- NULL
  phylo.d(compdata, binvar = trait, permut = 1)$DEstimate
}



tree_stats <- function(this_tree) {
  cat("\nAnalyzing: 0% [")
  if (any(is.na(this_tree))) {
    cat("----------]")
    return(NA)
  } else {
    
    
    
    ##### (0) Pull necessary variables from simulated trees and organize into a single object for all the tests below to pull from.
    
    #str(all_trees)
    #str(this_tree)
    
    
    ## 0a) Branch lengths
    Branch_Lengths <- this_tree$edge.length
    number_of_branches <- length(Branch_Lengths)
    
    # Anchor test = PD (Faith's phylogenetic diversity)
    Pylo_diversity_is_sum_of_BL <- sum(Branch_Lengths)
    
    # avPD -- Average phylogenetic diversity
    average_phylogenetic_diversity_is_mean_of_BL <- mean(Branch_Lengths)
    
    variance_Pylo_diversity_is_variance_of_BL <- var(Branch_Lengths)
    cat("-")
    ## 0b) Pairwise distance between tips
    Pairwise_dist <- cophenetic.phylo(this_tree)
    cat("-")
    # 2b) Pairwise distance -- Sum of pairwise distances
    
    # F -- Extensive quadratic entropy
    F_quadratic_entropy_is_sum_of_PD <- sum(Pairwise_dist)
    
    #Mean inter-species distances
    
    # Anchor test = MPD (mean pairwise distance)
    
    Mean_pairwise_distance <- mean(Pairwise_dist)
    
    cat("-")
    #Pairwise distance/all distances -- Variance of pairwise distances
    
    # Anchor test = VPD (variation of pairwise distance)
    
    variance_pairwise_distance <- var(as.vector(Pairwise_dist))
    
    
    
    
    ## 0c) Phylogenetic isolation
    
    # Using equal.splits method, faster computation
    Evolutionary_distinctiveness <- evol.distinct2(this_tree, type = "fair.proportion")
    
    cat("-")
    # ED - Summed evolutionary distinctiveness
    
    Evolutionary_distinctiveness_sum <- sum(Evolutionary_distinctiveness)
    
    ## 3d) Phylogenetic isolation -- Mean of species evolutionary distinctiveness
    
    # mean(ED)
    
    mean_Phylogenetic_isolation <- mean(Evolutionary_distinctiveness)
    
    ## 4d) Phylogenetic isolation -- Variance of species isolation metrics
    
    #var(ED)
    
    variance_Phylogenetic_isolation <- var(Evolutionary_distinctiveness)
    cat("-")
    
    ## Tree topology
    
    #Gamma index
    
    ltts <- ltt(this_tree, gamma = TRUE, plot = FALSE)
    lineages_through_time <- as.numeric(ltts[[1]])
    time_steps <- as.numeric(ltts[[2]])
    gamma <- ltts[[3]]
    gamma_p_value <- ltts[[4]]
    cat("-")
    
    colless_stat <- colless(as.treeshape(this_tree))
    sackin_index <- sackin(as.treeshape(this_tree))
    tree_shape_stat <- shape.statistic(as.treeshape(this_tree))
    
    ##### (5) Tree metric -- Macroevolutionary - Rate and rate changes ###############
    ##################################################
    
    ## Speciation vs extinction rates and Net diversification
    bds <- bd(this_tree)
    speciation_rate <- bds[1]
    extinction_rate <- bds[2]
    extinction_per_speciation <- bds[3]
    speciation_minus_extinction <- bds[4]
    cat("-")
    
    
    
    
    
    results_summary_matrix_1 <- cbind(
      
      number_of_branches,
      Pylo_diversity_is_sum_of_BL,
      average_phylogenetic_diversity_is_mean_of_BL,
      variance_Pylo_diversity_is_variance_of_BL,
      
      F_quadratic_entropy_is_sum_of_PD,
      Mean_pairwise_distance,
      variance_pairwise_distance,
      
      colless_stat ,
      sackin_index ,
      tree_shape_stat,
      
      Evolutionary_distinctiveness_sum,
      mean_Phylogenetic_isolation,
      variance_Phylogenetic_isolation,
      
      gamma,
      gamma_p_value,
      speciation_rate,
      extinction_rate,
      extinction_per_speciation,
      speciation_minus_extinction
    )
    rownames(results_summary_matrix_1) <- 1
    
    results_summary_matrix_2 <- cbind(
      c(Evolutionary_distinctiveness,NA),
      lineages_through_time,
      time_steps
    )
    colnames(results_summary_matrix_2) <- c("Evolutionary_distinctiveness",
                                            "lineages_through_time", "time_steps")
    head(results_summary_matrix_2)
    
    ### Returns from function in list form
    returns <- list(
      #Branch_Lengths,
      #Pairwise_dist,
      results_summary_matrix_1,
      results_summary_matrix_2
      
    )
    
    names(returns) <- c(
      #"Branch_Lengths",
      #"Pairwise_distance",
      "results_summary_of_single_value_outputs",
      "results_summary_matrix_of_multi_value_outputs"
    )
    cat("] 100%")
    
    return(returns)
    
  }
}

####### Bruno #######

run_sampler <- function (x, n, alpha, n_start = 1, return_start = FALSE,
                         starting = NULL) {
  
  # Error control
  .error_control(x, n, alpha, n_start, starting, return_start)
  x <- (x - min(x)) / (max(x) - min(x))
  # Names
  row_names_x <- rownames(x)
  
  # Starter
  n_row <- nrow(x)
  if (any(is.null(starting))) {
    starter <- sample(1:n_row, n_start)
    first <- row_names_x[starter]
  } else {
    n_start <- length(starting)
    first <- starting
    starter <- which(row_names_x %in% starting)
  }
  selected <- character(n)
  
  if (n_start == 1) {
    # Sample second
    dist_rela <- x[starter, -starter] ^ alpha
    if (any(is.infinite(dist_rela))) {
      stop("alpha is too high, infinite values generated.")
    }
    prob <- dist_rela/sum(dist_rela)
    second <- sample(names(prob), 1, prob = prob)
    selected[1:2] <- c(first, second)
    start_loop <- 3
  } else {
    start_loop <- n_start + 1
    selected[1:n_start] <- first
  }
  # loop
  for (i in start_loop:n) {
    positions <- row_names_x %in% selected
    dist_rela <- apply(x[positions, !positions], 2, min) ^ alpha
    if (any(is.infinite(dist_rela))) {
      stop("alpha is too high, infinite values generated.")
    }
    prob1 <- dist_rela / sum(dist_rela)
    dist_rela2 <- apply(x[positions, !positions], 2, min) ^ alpha
    if (any(is.infinite(dist_rela2))) {
      stop("alpha is too high, infinite values generated.")
    }
    prob2 <- dist_rela2 / sum(dist_rela2)
    prob <- (prob2 + prob1) / 2
    selected[i] <- sample(names(prob), 1, prob = prob)
  }
  if (return_start) {
    return(list('Sampling_selection' = selected, "Starting_points" = first))
  } else {
    return(selected)
  }
}



# Error control function
.error_control <- function (x, n, alpha, n_start, starting, return_start) {
  if (!is.matrix(x)) {
    stop("x must be a matrix.")
  }
  if (nrow(x) != ncol(x)) {
    stop("x must be a symmetric distance matrix.")
  }
  if (is.null(colnames(x)) | is.null(rownames(x))) {
    stop("x should have both columns and rows named.")
  }
  if (any(duplicated(rownames(x)))) {
    stop("duplicated names in x.")
  }
  if (!is.numeric(n)) {
    stop("n must be a number")
  }
  if (length(n) != 1) {
    stop("n should have length 1")
  }
  if (n < 2) {
    stop("n must be higher than 1")
  }
  if (n/n != 1) {
    stop("n must be an integer")
  }
  if (!is.numeric(alpha)) {
    stop("alpha should a number")
  }
  if (length(alpha) != 1) {
    stop("alpha should have length 1")
  }
  if (!is.numeric(n_start)) {
    stop("n_start must be a number")
  }
  if (length(n_start) != 1) {
    stop("n_start should have length 1")
  }
  if (n_start < 1) {
    stop("n_start must be a positive value equals or bigger than 1")
  }
  if (n_start/n_start != 1) {
    stop("n_start must be an integer")
  }
  if (!is.logical(return_start)){
    stop("return_start must be logical")
  }
  invisible(NULL)
}


run_sampler_phy <- function (x, n, alpha, dist.func = cophenetic,
                             n_start = 1, return_start = FALSE,
                             starting = NULL) {
  #
  if (class(x) != 'phylo') {
    stop("x must be a phylogeny")
  }
  if (class(dist.func) != "function") {
    stop("dist.func must be a function")
  }
  tree <- x
  x <- dist.func(tree)
  if(!all(tree$tip.label %in% rownames(x))) {
    stop("dist.func provided do not result in a distance matrix
         with names corresponding to the tree tip labels.")
  }
  selected <- run_sampler(x, n, alpha, n_start, return_start,
                          starting)
  if (return_start) {
    selected[[1]] <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% selected[[1]]])
    return(selected)
  } else {
    selected <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% selected])
    return(selected)
  }
  }




####################
#### file pulls

avail_files <- list.files(path="~/Desktop/Finished_tree_sets", recursive=TRUE, pattern="Full_trees*")

finished_files <- list.files(path="~/Desktop/Finished_tree_prunes", recursive=TRUE)

for(i in 1:length(finished_files)){
  trunker <- strsplit(finished_files[[i]], "_")[[1]]
  finished_files[i] <- paste(trunker[-1:-2], collapse ="_")
}
#finished_files

already_done <- which(!(avail_files %in% finished_files))

if(length(already_done) != 0)(unique_remaining <- avail_files[already_done])
#length(unique_remaining)

args <- commandArgs(trailingOnly = FALSE)
cluster_value <- as.numeric(args[7])

pre_files <-  sample(unique_remaining, 1)


files <- paste0("~/Desktop/Finished_tree_sets/", pre_files)
load(as.character(files))

### Prunes
percent_remove <- as.integer(runif(1, min=1, max=100))

n<-length(full_tree$tip.label)
ee<-setNames(full_tree$edge.length[sapply(1:n,function(x,y)   which(y==x),y=full_tree$edge[,2])],full_tree$tip.label)
trim <- length(ee) * (percent_remove/100)
dtips <- names(sample(ee, size= trim)) 
trimmer <- which(full_tree$tip.label %in% dtips)
random_trimmed_tree <-  drop.tip(full_tree, trimmer)
random_trimmed_tree


# proportionally shallower splits
n<-length(full_tree$tip.label)
ee<-setNames(full_tree$edge.length[sapply(1:n,function(x,y)   which(y==x),y=full_tree$edge[,2])],full_tree$tip.label)
trim <- length(ee) * (percent_remove/100)
dtips <- names(sample(ee, size= trim, prob=as.numeric(ee))) 
trimmer <- which(full_tree$tip.label %in% dtips)
shallow_trimmed_tree <-  drop.tip(full_tree, trimmer)
shallow_trimmed_tree

# Proportionally deeper splits
n<-length(full_tree$tip.label)
ee<-setNames(full_tree$edge.length[sapply(1:n,function(x,y)   which(y==x),y=full_tree$edge[,2])],full_tree$tip.label)
trim <- length(ee) * (percent_remove/100)
dtips <- names(sample(ee, size= trim, prob=1/as.numeric(ee)))
trimmer <- which(full_tree$tip.label %in% dtips)
deep_trimmed_tree <-  drop.tip(full_tree, trimmer)
deep_trimmed_tree

trim <- 1000 * (percent_remove/100)
dtips <- run_sampler_phy(full_tree, trim, alpha = 50)
trimmer <- which(full_tree$tip.label %in% dtips$tip.label)
overdispersed_tree <-  drop.tip(full_tree, trimmer)
overdispersed_tree

trim <- 1000 * (percent_remove/100)
dtips <- run_sampler_phy(full_tree, trim, alpha = -50)
trimmer <- which(full_tree$tip.label %in% dtips$tip.label)
underdispersed_tree <-  drop.tip(full_tree, trimmer)
underdispersed_tree

require(ape)
pruned_trees<-list(full_tree, random_trimmed_tree, shallow_trimmed_tree, deep_trimmed_tree, overdispersed_tree, underdispersed_tree)
class(pruned_trees)<-"multiPhylo"




stats_list_1 <- tree_stats(pruned_trees[[1]])
stats_list_2 <- tree_stats(pruned_trees[[2]])
stats_list_3 <- tree_stats(pruned_trees[[3]])
stats_list_4 <- tree_stats(pruned_trees[[4]])
stats_list_5 <- tree_stats(pruned_trees[[5]])
stats_list_6 <- tree_stats(pruned_trees[[6]])
Summary_stats_matrix <- matrix(rep(NA,19), 1, 19, byrow=TRUE)
Summary_stats_matrix <- rbind(stats_list_1$results_summary_of_single_value_outputs, stats_list_2$results_summary_of_single_value_outputs, stats_list_3$results_summary_of_single_value_outputs, stats_list_4$results_summary_of_single_value_outputs, stats_list_5$results_summary_of_single_value_outputs,
                              stats_list_6$results_summary_of_single_value_outputs)
Prune_type <- c("full_tree", "random_trimmed_tree", "shallow_trimmed_tree", "deep_trimmed_tree", "overdispersed_tree", "underdispersed_tree")
Summary_stats_matrix <- as.data.frame(cbind(Prune_type, Summary_stats_matrix))
rownames(Summary_stats_matrix) <- c(1:6)
#Summary_stats_matrix



BL_matrix <- matrix(rep(NA, 12000), 2000, 6)
BL_matrix[1:length(pruned_trees[[1]]$edge[,2]),1] <- pruned_trees[[1]]$edge[,2]
BL_matrix[1:length(pruned_trees[[2]]$edge[,2]),2] <- pruned_trees[[2]]$edge[,2]
BL_matrix[1:length(pruned_trees[[3]]$edge[,2]),3] <- pruned_trees[[3]]$edge[,2]
BL_matrix[1:length(pruned_trees[[4]]$edge[,2]),4] <- pruned_trees[[4]]$edge[,2]
BL_matrix[1:length(pruned_trees[[5]]$edge[,2]),5] <- pruned_trees[[5]]$edge[,2]
BL_matrix[1:length(pruned_trees[[6]]$edge[,2]),6] <- pruned_trees[[6]]$edge[,2]
#rownames(BL_matrix) <- pruned_trees[[1]]$edge[,1]
colnames(BL_matrix) <- Prune_type
#BL_matrix



time_stats_matrix <- matrix(rep(NA, 1001*18), 1001, 18)
time_stats_matrix[1:length(stats_list_1$results_summary_matrix_of_multi_value_outputs[,1]),1:3] <- stats_list_1$results_summary_matrix_of_multi_value_outputs
time_stats_matrix[1:length(stats_list_2$results_summary_matrix_of_multi_value_outputs[,1]),4:6] <- stats_list_2$results_summary_matrix_of_multi_value_outputs
time_stats_matrix[1:length(stats_list_3$results_summary_matrix_of_multi_value_outputs[,1]),7:9] <- stats_list_3$results_summary_matrix_of_multi_value_outputs
time_stats_matrix[1:length(stats_list_4$results_summary_matrix_of_multi_value_outputs[,1]),10:12] <- stats_list_4$results_summary_matrix_of_multi_value_outputs
time_stats_matrix[1:length(stats_list_5$results_summary_matrix_of_multi_value_outputs[,1]),13:15] <- stats_list_5$results_summary_matrix_of_multi_value_outputs
time_stats_matrix[1:length(stats_list_6$results_summary_matrix_of_multi_value_outputs[,1]),16:18] <- stats_list_6$results_summary_matrix_of_multi_value_outputs

Prune_type_names <- c("full_tree_Evolutionary_distinctiveness","full_tree_lineages_through_time","full_tree", "random_trimmed_tree_Evolutionary_distinctiveness","random_trimmed_tree_lineages_through_time","random_trimmed_tree_time_steps", "shallow_trimmed_tree_Evolutionary_distinctiveness","shallow_trimmed_tree_lineages_through_time","shallow_trimmed_tree_time_steps", "deep_trimmed_tree_Evolutionary_distinctiveness","deep_trimmed_tree_lineages_through_time","deep_trimmed_tree_time_steps", 
                      "overdispersed_tree_Evolutionary_distinctiveness", "overdispersed_tree_lineages_through_time","overdispersed_tree_time_steps",
                      "underdispersed_tree_Evolutionary_distinctiveness", "underdispersed_tree_lineages_through_time", "underdispersed_tree_time_steps")
colnames(time_stats_matrix) <- Prune_type_names
#time_stats_matrix


pruned_trees_and_stats <- list(paste0("Percent removed ", percent_remove), full_tree, random_trimmed_tree, shallow_trimmed_tree, deep_trimmed_tree, overdispersed_tree, underdispersed_tree, Summary_stats_matrix, time_stats_matrix, BL_matrix)
pruned_trees_and_stats


save(pruned_trees_and_stats, file=paste("~/Desktop/Finished_tree_prunes/PercentRemoved_", percent_remove, "_", pre_files, sep=""))


