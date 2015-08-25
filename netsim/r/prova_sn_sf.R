setwd("G:\\R\\r2c\\R2Clib")

source("simulatenet.R")

ris <- simulatenet_dll(N=5, connectivity="scale free", times=seq(1,5))
cat("ok")