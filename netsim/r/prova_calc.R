#setwd("D:\\Progetti R\\netsim_hg\\sim")
setwd("/home/marco/Documenti/Progetti_R/netsim/sim")

#dyn.load("calc.dll")
dyn.load("calc.so")

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

t <- seq(0, 10, 0.1)
y <- calcv_dll("pulse(x, 0.5, 10)", t)
y
plot(t, y, type="l")

t <- 1
y1 <- calc_dll("if(15.5 > 0, 3, 0)", t)
y1
