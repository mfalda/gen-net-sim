setwd("G:\\R\\r2c\\R2Clib")

dyn.load("mandel.dll")

mandel_dll <- function(param, n)
{
    stopifnot(is.vector(param))
    stopifnot(is.numeric(n))
    matr <- .Call("mandel",
                 as.vector(param),
                 as.integer(n)
					)
    return(matr)
}

min_x <- 0.32
max_x <- 0.38
min_y <- -0.38
max_y <- -0.32
risol <- 50
julia <- 0 + 0i
n <- 500
param <- c(min_x, max_x, min_y, max_y, risol, Re(julia), Im(julia))

cat("Calcolo la matrice ... ")
t <- system.time(m1 <- mandel_dll(param, n))
cat(sprintf("fatto: (%3.3f secondi)\n", t[1][[1]]))
passo_x <- (max_x - min_x) / risol
passo_y <- (max_y - min_y) / risol
image(seq(min_x, max_x, passo_x), seq(min_y, max_y, passo_y), m1, xlab="Real", ylab="Imag", col=c("black", rainbow(n)))
