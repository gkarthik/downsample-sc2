#!/usr/bin/env R

library(ape)

## 29891 columns in aligment
args = commandArgs(trailingOnly=TRUE)

tree_file <- args[1]
out_file <- args[2]

min_length <- 0.99999/29891
tree <- read.tree(tree_file)

collapsed_tree <- di2multi(tree, tol = min_length)
collapsed_tree <- root(collapsed_tree, outgroup = "EPI_ISL_402124|2019-12-30|China|Hubei|Wuhan", resolve.root = TRUE)
write.tree(out_file)
