
dyn.load("createRules.so")

create.Rules <-
    function (matr)
{
    stopifnot(is.matrix(matr))
    mat <- .Call("createRules",
                 as.matrix(matr))
    return(mat)
}

m <- matrix(seq(10,15), 2, 3)
m
create.Rules(m)
