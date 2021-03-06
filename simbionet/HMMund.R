HMMund<-function(N=50, Cf.cl=0.3, gamma=2.2, DEGREE=NULL, MODULES, prior.p.subnet=NULL, max.con=12, sepgraph=TRUE, r.tol=0.1, a.tol=1, iter=1)
{
    if (!is.numeric(N))
        stop("`N' must be numeric")

    if (!is.numeric(Cf.cl))
        stop("`max.reg' must be numeric")

    if (!is.numeric(gamma))
        stop("`gamma' must be numeric")

    if (!is.null(DEGREE)) {
        if (!is.numeric(DEGREE))
            stop("DEGREE must be numeric")
        if (length(which(DEGREE<0))>0)
            stop("DEGREE must contain positive numbers")
        if (sum(DEGREE,na.rm=TRUE)!=N)
            stop("sum(DEGREE) must be equal to N")
        if (length(DEGREE)!=(N+1))
            stop("DEGREE must be of length (N+1)")
        if (max(DEGREE,na.rm=TRUE)>N)
            stop("DEGREE must have AT MAXIMUM value N")
    }

    if (!is.list(MODULES))
        stop("'MODULES' must be a list")

    L<-length(MODULES)

    if (is.null(prior.p.subnet))
        prior.p.subnet<-rep(1/L,L)
    if (!is.numeric(prior.p.subnet))
        stop("`prior.p.subnet' must be numeric")
    if (length(prior.p.subnet)!=L)
        stop("`prior.p.subnet' must have length L=",L)
    if (length(which(is.na(prior.p.subnet)))>0)
        stop("`prior.p.subnet' can't contain NA values")
    if ((length(which(prior.p.subnet>1))>0)|(length(which(prior.p.subnet<0))>0))
        stop("`prior.p.subnet' must contain values greater than 0 and lower than 1")
    if (round((sum(prior.p.subnet))-1,10)!=0)
        stop("`prior.p.subnet' must have sum 1")


    if (!is.numeric(max.con))
        stop("`max.con' must be numeric")


	if (!is.logical(sepgraph))
        stop("'sepgraph' must be TRUE or FALSE")

    if (!is.numeric(r.tol))
        stop("`r.tol' must be numeric")
    else {
        if ((r.tol>1)|(r.tol<0))
            stop("`r.tol' must be greater than 0 and lower than 1")
    }

    if (!is.numeric(a.tol))
        stop("`a.tol' must be numeric")
    else {
        if (a.tol<0)
            stop("`a.tol' must be a positive number")
    }

    ris <- .Call("HMMund",
        as.integer(N),
        as.numeric(Cf.cl),
        as.numeric(gamma),
        as.vector(DEGREE),
        as.list(MODULES),
        as.vector(prior.p.subnet),
        as.integer(max.con),
        as.logical(sepgraph),
        as.numeric(r.tol),
        as.numeric(a.tol),
        as.integer(iter)
    )

    return (ris)
}