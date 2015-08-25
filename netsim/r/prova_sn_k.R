#setwd("D:\\Progetti R\\netsim_hg\\sim")
setwd("/home/marco/Documenti/Progetti_R/netsim/sim")

source("simulatenet.R")

net<-simulatenet(N=20,connectivity="geometric",kappa=5, method="Euler")
#net <- simulatenet(N=50, ,connectivity="MTM", f.pr.and="1 / (1 + _e ^ (-10 * (x - 0.5)))", act.fun="sigmoidal", method="rkf45")
cat("ok")