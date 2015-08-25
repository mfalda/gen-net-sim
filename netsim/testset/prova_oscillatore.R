setwd("G:\\R\\r2c\\R2Clib")

library(odesolve)

dyn.load("lsoda_oscillatore.dll")

# oscillatore: m * x'' + b * x' + k * x = F cos(o * t)
# derivate: [y[1], -(k / m) * y[0] - (b / m) * y[1] + F * cos(o * t)]
f <- function(t, y, parms)  
{
	# m, b, k, F, omega
	m <- parms[[1]]	
	b <- parms[[2]]
	k <- parms[[3]]
	F <- parms[[4]]
	omega <- parms[[5]]
	dydt<- rbind(y[2], -(k / m) * y[1] - (b / m) * y[2] + F * cos(omega * t));
	return (list(dydt))
}

oscillatore_dll <- function(args, x0, times, method="rkf45", atol=1e-6, rtol=0.0)
{
    stopifnot(is.list(args))
    stopifnot(is.vector(x0))
    stopifnot(is.vector(times))
    stopifnot(is.character(method))
    stopifnot(is.numeric(atol))
    stopifnot(is.numeric(rtol))
    matr <- .Call("lsoda_oscillatore",
                 as.list(args),
                 as.vector(x0),
                 as.vector(times),
                 as.character(method),
                 as.numeric(atol),
                 as.numeric(rtol)
					)
    return(matr)
}

# m, b, k , F, omega
m = 5
b = 2
k = 5
F = 10
omega = 3.1415
p <- list(m, b, k, F, omega)
x0 <- c(1, 0.0)
t <- seq(0, 10, 0.1)
ris <- oscillatore_dll(p, x0, t, method="rkck")
ris1 <- lsoda(x0, t, f, p, rtol=1e-10, atol=1e-10)
z <- sqrt(b ^ 2 + (omega * m - k / m) ^ 2)
phi <- atan((omega * m - k / m) / b)
x <- function(t) {return (F / z / omega * sin(omega * t - phi))}
plot(ris[,1],abs(ris[,2] - ris1[,2]), type='l')