setwd("/home/marco/R2Clib")

dyn.load("lsoda_oscillatore.so")

oscillatore_dll <- function(args, x0, times)
{
    stopifnot(is.list(args))
    stopifnot(is.vector(x0))
    stopifnot(is.vector(times))
    matr <- .Call("lsoda_oscillatore",
                 as.list(args),
                 as.vector(x0),
                 as.vector(times)
					)
    return(matr)
}

# m, b, k , F, omega
m = 5
b = 0.5
k = 5
F = 0
omega = 3.1415 * 2
p <- list(m, b, k, F, omega)
x0 <- c(-0.05, 0.0)
t <- seq(0, 100, 1)
ris <-oscillatore_dll(p, x0, t)
z <- sqrt(b ^ 2 + (omega * m - k / m) ^ 2)
phi <- atan((omega * m - k / m) / b)
x <- function(t) {return (F / z / omega * sin(omega * t - phi))} 
plot(ris[,1],ris[,2], type='l')