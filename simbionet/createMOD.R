createMOD<-function(m=4, auto=FALSE)
{
    if (!is.numeric(m))
        stop("'m' must be numeric")

    if ((m<2) | (m>5))
        stop("'m' must be greater than 1 and lower than 6")

    if (!is.logical(auto))
        stop("'auto' must be TRUE or FALSE")

    ris <- .Call("createMOD",
        as.integer(m),
        as.integer(auto)
    )

    return (ris)
}
