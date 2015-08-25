setwd("/Users/mf/Desktop/linux/sim")

dyn.load("calc.so")

calc_dll <- function(nome, args)
{
    stopifnot(is.character(nome))
    stopifnot(is.vector(args))
    vett <- .Call("calc",
                 as.character(nome),
                 as.character(args)
					)
    return(vett)
}

calcv_dll <- function(nome, args)
{
    stopifnot(is.character(nome))
    stopifnot(is.vector(args))
    vett <- .Call("calcv",
                 as.character(nome),
                 as.character(args)
					)
    return(vett)
}

t <- seq(0, 1, 0.1)
t
cat(calcv_dll("sin(x)", t), "\n")