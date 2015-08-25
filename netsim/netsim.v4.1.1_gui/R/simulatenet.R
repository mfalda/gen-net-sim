simulatenet<-function(N=50,connectivity="MTM",max.reg=12,gamma=2.2,INdegree="free", Cf.cl=0.4,num.subnet=c(5,5,10),kappa=3, f.pr.and="", Xmin=NULL, Xmax=NULL,lambda=c(1,0.1),x0=c(0,1),weight.par=c(1,0),act.fun="sigmoidal", alpha=c(10,0.2), beta=c(0.5,0.01),times=seq(0,5,0.05),method="rkf45",params=c(2, 2,2,0,0,1),save=FALSE,itera=1)
{
    if (!is.numeric(N))
        stop("`N' must be numeric")

    if (!is.numeric(max.reg))
        stop("`max.reg' must be numeric")

    if (!is.numeric(gamma))
        stop("`gamma' must be numeric")

    if (!is.numeric(Cf.cl))
        stop("`max.reg' must be numeric")

    if (!is.numeric(num.subnet))
        stop("`num.subnet' must be numeric")

    if (!is.numeric(kappa))
        stop("`kappa' must be numeric")

    if (!is.null(f.pr.and)) {
        if (!is.character(f.pr.and))
            stop("`f.pr.and' must be NULL or a string")
    }
    else
        f.pr.and <- ""

    if (is.null(Xmin))
        Xmin<-rep(0,N)
    else if (!is.numeric(Xmin))
        stop("`Xmin' must be numeric")

    if (is.null(Xmax))
        Xmax<-rep(10,N)
    else if (!is.numeric(Xmax))
        stop("`Xmax' must be numeric")

    if (!is.numeric(lambda))
        stop("`lambda' must be numeric")

    if (!is.numeric(x0))
        stop("`x0' must be numeric")

    if (!is.numeric(weight.par))
        stop("`weight.mean' must be numeric")

    if (!is.numeric(alpha))
        stop("`alpha' must be numeric")

    if (!is.numeric(beta))
        stop("`beta' must be numeric")

    if (!is.numeric(times))
        stop("`times' must be numeric")

    if (!is.numeric(params))
        stop("`params' must be numeric")

    if (!is.logical(save))
        stop("`save' must be logical")

    if (!is.numeric(itera))
        stop("`itera' must be numeric")

    ris <- .Call("simulatenet",
        as.integer(N),
        as.character(connectivity),
        as.integer(max.reg),
        as.numeric(gamma),
        as.character(INdegree),
        as.numeric(Cf.cl),
        as.vector(num.subnet),
        as.numeric(kappa),
        as.character(f.pr.and),
        as.vector(Xmin),
        as.vector(Xmax),
        as.vector(lambda),
        as.vector(x0),
        as.vector(weight.par),
        as.character(act.fun),
        as.vector(alpha),
        as.vector(beta),
        as.vector(times),
        as.character(method),
        as.vector(params),
        as.logical(save),
        as.integer(itera)
    )
    return(ris)
}
