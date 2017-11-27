library(diversitree)
#setwd("~/Desktop")


# make tree
params <- c(1, 1, 1, 1, # lambda 00, 10, 01, 11
0.3, 0.3, 0.3, 0.3,  # mu 00, 10, 01, 11
0.14, 0.12, 0, #q00.10, q00.01, q00.11
0.08, 0, 0.05, #q10.00, q10.01, q10.11
0.11, 0, 0.07, #q01.00, q01.10, q01.11
0, 0.13, 0.06) #q11.00, q11.10, q11.01
names <- c("lambda1",    "lambda2",    "lambda3" ,   "lambda4",        "mu1" ,       "mu2"  ,      "mu3",        "mu4"  ,      "q12" ,       "q13",        "q14",        "q21" ,       "q23"    ,    "q24",       "q31"  ,      "q32"  ,  "q34"  ,      "q41"  ,      "q42"      ,  "q43" )

origSize <- 1000
full_tree <- tree.musse(params, max.taxa=origSize, include.extinct =FALSE, x0=4)

args <- commandArgs(trailingOnly = FALSE)
NAI <- as.numeric(args[7])
replicate <- as.integer(runif(1, 1, 100000))
param_file_name <- paste(paste(names, params,  sep="_" ), collapse="_")
file_name <- paste("./Full_trees/Full_tree","tips", origSize , param_file_name,"rep", replicate, ".Rdata", sep="_")

save(full_tree, file= file_name)
