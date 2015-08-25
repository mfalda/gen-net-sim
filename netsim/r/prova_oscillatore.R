setwd("D:\\Progetti R\\netsim_hg\\sim")

dyn.load("lsoda_oscillatore.dll")
#dyn.unload("lsoda_oscillatore.dll")

oscillatore_dll <- function(args, x0, times, method="rkf45", atol=1e-6, rtol=0.0, stat_t=0.001, stat_w=1)
{
    stopifnot(is.list(args))
    stopifnot(is.vector(x0))
    stopifnot(is.vector(times))
    stopifnot(is.character(method))
    stopifnot(is.numeric(atol))
    stopifnot(is.numeric(rtol))
    stopifnot(is.numeric(stat_w))
    stopifnot(is.numeric(stat_t))
    matr <- .Call("lsoda_oscillatore",
                 args,
                 x0,
                 times,
                 method,
                 atol,
                 rtol,
                 stat_t,
                 stat_w
					)
    return(matr)
}

# m, b, k , F, omega
m = 5
b = 1
k = 5
F = 0
omega = 3.1415
p <- list(m, b, k, F, omega)
x0 <- c(-0.05, 0.0)
t <- seq(0, 100, 1)
ris <-oscillatore_dll(p, x0, t, method="rkck", stat_t=0.1, stat_w=0.3)
z <- sqrt(b ^ 2 + (omega * m - k / m) ^ 2)
phi <- atan((omega * m - k / m) / b)
x <- function(t) {return (F / z / omega * sin(omega * t - phi))}
plot(ris[,1],ris[,2], type='l')