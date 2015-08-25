test <- "prova_sn_actfun"
ris <- as.matrix(read.table(sprintf("G:\\R\\r2c\\R2Clib\\testset\\%s\\R\\SIMdata1.txt", test)))
ris1 <- as.matrix(read.table(sprintf("G:\\R\\r2c\\R2Clib\\testset\\%s\\distrib\\SIMdata1.txt", test)))
mx <- 0
for (i in 1:5) {
  diff <- abs(ris1[i,] - ris[i,])
  m1 <- max(diff)
  cat("Profilo", i, " = ", m1, "\n")
  if (m1 > mx)
    mx <- m1
}
cat(sprintf("La massima differenza è %.e\n", mx))
m2 <- max(abs(ris[,1] - ris1[,1]))
if (m2 > 0)
  cat(sprintf("Gli X0 differiscono al massimo di %.e\n", m2))
plot(1:51, ris[1,])
lines(1:51, ris1[1,], type="l", col="red")