source("simulatenet.R")

system.time(simulatenet(N=5, connectivity="MTM", method="rkf45", itera=3))
cat("ok")