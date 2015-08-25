setwd("D:\\Progetti R\\netsim_hg\\sim")

dyn.load("mandel_calc.dll")

mandelc_dll <- function(param, reale, img, n)
{
    stopifnot(is.vector(param))
    stopifnot(is.character(reale))
    stopifnot(is.character(img))
    stopifnot(is.numeric(n))
    matr <- .Call("mandel",
                 param,
                 reale,
                 img,
                 n
					)
    return(matr)
}

min_x <- 0.32
max_x <- 0.38
min_y <- -0.38
max_y <- -0.32
risol <- 90
julia <- 0 + 0i
n <- 5000
param <- c(min_x, max_x, min_y, max_y, risol, Re(julia), Im(julia))

cat("Calcolo la matrice ... ")
t <- system.time(m1 <- mandelc_dll(param, "x * x - y * y + z", "2 * x * y + z", n))
cat(sprintf("\nfatto: (%3.3f secondi)\n", t[1][[1]]))
passo_x <- (max_x - min_x) / risol
passo_y <- (max_y - min_y) / risol
image(seq(min_x, max_x, passo_x), seq(min_y, max_y, passo_y), m1, xlab="Real", ylab="Imag", col=c("black", rainbow(n)))
