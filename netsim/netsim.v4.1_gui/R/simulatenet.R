simulatenet <- function(N=50, connectivity="MTM", max.reg=12, gamma=2.2, INdegree="free", Cf.cl=0.4, num.subnet=c(5,5,10), kappa=3, f.pr.and="", act.fun="sigmoidal", alpha=c(10,0.2), beta=c(0.5,0.01), lambda=c(1,0.1), Xmin=NULL, Xmax=NULL, X0=c(0,1), weight.par=c(1,0), params=c(2,2,2,0,0,1), times=seq(0,5,0.05), stat_thr=0.001, stat_width=0, method="rkf45", itera=1, save=FALSE)
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

    if (!is.numeric(X0))
        stop("`X0' must be numeric")

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

    if (!is.numeric(stat_thr))
        stop("'stat_thr' must be numeric")

    if (!is.numeric(stat_width))
        stop("'stat_width' must be numeric")

    ris <- .Call("simulatenet",
        N,
        connectivity,
        max.reg,
        gamma,
        INdegree,
        Cf.cl,
        num.subnet,
        kappa,
        f.pr.and,
        act.fun,
        alpha,
        beta,
        lambda,
        Xmin,
        Xmax,
        X0,
        weight.par,
        params,
        times,
        stat_thr,
        stat_width,
        method,
        itera,
        save
    )
    return(ris)
}
