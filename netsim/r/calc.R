setwd("D:\\Progetti R\\netsim_hg\\sim")

dyn.load("calc.dll")

calc_dll <- function(nome, args)
{
    stopifnot(is.character(nome))
    stopifnot(is.vector(args))
    vett <- .Call("calc",
                 nome,
                 args
					)
    return(vett)
}

calcv_dll <- function(nome, args)
{
    stopifnot(is.character(nome))
    stopifnot(is.vector(args))
    vett <- .Call("calcv",
                 nome,
                 args
					)
    return(vett)
}

t <- seq(0, 10, 0.01)
t
y <- calcv_dll("periodic1(x, 1, 0.15, 0.95)", t)
y
plot(t, y, type="l")

t <- 1
y1 <- calc_dll("15.5 % 3", t)
y1