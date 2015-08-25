setwd("/home/marco/R2Clib")

dyn.load("simulatenet.so")

simulatenet_dll<-function(N=50,connectivity="MTM",max.reg=12,gamma=2.2,INdegree="free", Cf.cl=0.4,num.subnet=c(5,5,10),kappa=3, f.pr.and="", Xmax=NULL,lambda=NULL,x0=NULL,weight.par=c(1,0),act.fun="sigmoidal", alpha=NULL, beta=NULL,times=seq(0,5,0.1),method="Euler",save=FALSE,ind.itera=1)
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
            stop("`f.pr.and' must be NULL or a function")
    }
    else
        f.pr.and <- ""

    if (is.null(Xmax))
        Xmax<-rep(10,N)
    else if (!is.numeric(Xmax))
        stop("`Xmax' must be numeric")

    if (is.null(lambda))
        lambda<-abs(rnorm(N,1,0.1))
    else if (!is.numeric(lambda))
        stop("`lambda' must be numeric")

    if (is.null(x0))
        x0<-runif(N)
    else if (!is.numeric(x0))
        stop("`x0' must be numeric")

    if (!is.numeric(weight.par))
        stop("`weight.mean' must be numeric")

    if (is.null(alpha))
        alpha<-abs(rnorm(N,10,0.2))
    else if (!is.numeric(alpha))
        stop("`alpha' must be numeric")

    if (is.null(beta))
        beta<-abs(rnorm(N,0.5,0.01))
    else if (!is.numeric(beta))
        stop("`beta' must be numeric")

    if (!is.numeric(times))
        stop("`times' must be numeric")

    if (!is.logical(save))
        stop("`save' must be logical")

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
        as.vector(Xmax),
        as.vector(lambda),
        as.vector(x0),
        as.vector(weight.par),
        as.character(act.fun),
        as.vector(alpha),
        as.vector(beta),
        as.vector(times),
        as.character(method),
        as.logical(save),
        as.integer(ind.itera)
    )
    return(ris)
}

ris <- simulatenet_dll(350,times=seq(1,5))
cat("ok")