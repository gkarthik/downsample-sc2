#!/usr/bin/env Rscript

library(ape)

## 29891 columns in aligment
args = commandArgs(trailingOnly=TRUE)

tree_file <- args[1]
out_file <- args[2]

tree <- read.tree(tree_file)

drop.tip(tree, "EPI_ISL_402124|2019-12-30|China|Hubei|Wuhan")

write.tree(tree, out_file)
