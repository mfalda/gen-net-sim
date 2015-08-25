setwd("D:\\Progetti R\\netsim_hg\\sim")

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
risol <- 500
julia <- 0 + 0i
param <- c(min_x, max_x, min_y, max_y, risol, Re(julia), Im(julia))
tempi <- c()
cat("Calcolo la matrice ... \n")
for (n in seq(100, 1000, 100)) { 
	for (i in seq(1, 10)) { 
		t <- system.time(m1 <- mandel_dll(param, n))
		#cat(sprintf("%d\t%3.3f\n", i, t[1][[1]]))
		tempi <- c(tempi, t[1][[1]])
	}
	cat(sprintf("%3.3f\t%3.3f\n", mean(tempi), sd(tempi)))
}

passo_x <- (max_x - min_x) / risol
passo_y <- (max_y - min_y) / risol
image(seq(min_x, max_x, passo_x), seq(min_y, max_y, passo_y), m1, xlab="Real", ylab="Imag", col=c("black", rainbow(n)))
