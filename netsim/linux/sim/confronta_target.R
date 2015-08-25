setwd("H:\\R\\r2c\\ok")


library("netsim")
dyn.load("target.dll")

target_C <-
    function (n, r, m, N, ext_in, ext_fun, t)
{
    stopifnot(is.vector(n))
    stopifnot(is.list(r))
    stopifnot(is.matrix(m))
    stopifnot(is.matrix(N))
    stopifnot(is.matrix(ext_in))
    stopifnot(is.vector(ext_fun))
    stopifnot(is.numeric(t))
    vett <- .Call("target",
                 as.matrix(n),
                 as.list(r),
                 as.matrix(m),
                 as.matrix(N),
                 as.matrix(ext_in),
                 as.vector(ext_fun),
                 as.numeric(t)
                 )
    return(vett)
}

dm = 3 #round(runif(1) * 10) + 1
y = rep(0, dm)
for (i in 1:dm)
    y[i] = runif(1)
y

R = list()
R[[1]] = dm
R[[2]] = c(-2, 1, 3)
R[[3]] = 1

R

m = matrix(0, dm, dm)
for (i in 1:dm) {
    for (j in 1:dm)
        m[i, j] = runif(1)
}

N = matrix(0, dm, dm)
for (i in 1:dm) {
    for (j in 1:dm)
        N[i, j] = runif(1)
}

ext_in = matrix()

ext_fun = vector()

t = runif(1)
t

r1 = target_C(y, R, m ,N, ext_in, ext_fun, t)
r1

r = target(y, R, m ,N, ext_in, ext_fun, t)
r

cat("Elementi diversi: ", which(r != r1), "\n")

