setwd("/home/marco/R2Clib")

dyn.load("simulateprofiles.so")

simulateprofiles_dll<-function(weights, Rules=NULL, f.pr.and=NULL, x0=NULL, Xmax=NULL, lambda=NULL, act.fun="sigmoidal", alpha=NULL, beta=NULL, times=seq(0,5,0.1), method="Euler", EXT.IN=NA, EXT.FUN=list())
{
    if (!is.numeric(weights))
        stop("'weights' must be numeric")
    if (!is.matrix(weights))
        stop("'weights' must be a matrix")
	N<-dim(weights)[1]

    if (!is.null(Rules)) {
        if (!is.list(Rules))
            stop("'rules' must be a list or NULL")
        if (!is.numeric(Rules[[1]]))
            stop("'Rules' must contain numeric values")
    }

   if (!is.null(f.pr.and)) {
        if (!is.character(f.pr.and))
            stop("`f.pr.and' must be a string")
   }
   else
        f.pr.and <- ""

    if (is.null(lambda))
        lambda<-abs(rnorm(N,1,0.1))
    if (!is.numeric(lambda))
        stop("`lambda' must be numeric")

    if (is.null(x0))
        x0<-runif(N)
    if (!is.numeric(x0))
        stop("'x0' must be numeric")
    if ((max(x0)>1)|(min(x0)<0))
        stop("'x0' must be in the range [0,1]")

    if (is.null(Xmax))
        Xmax<-rep(10,N)
    if (!is.numeric(Xmax))
        stop("'Xmax' must be numeric")

    if (!is.numeric(lambda))
        stop("'lambda' must be numeric")

    if (is.null(alpha))
        alpha<-abs(rnorm(N,10,0.2))
    if (!is.numeric(alpha))
        stop("`alpha' must be numeric")

    if (is.null(beta))
        beta<-abs(rnorm(N,0.5,0.01))
    if (!is.numeric(beta))
        stop("`beta' must be numeric")

    if (!is.numeric(times))
        stop("'times' must be numeric")

    if (length(EXT.FUN)>0) {
        if (!is.numeric(EXT.IN))
            stop("'EXT.IN' must be numeric")
        if (!is.matrix(EXT.IN))
            stop("'EXT.IN' must be a matrix")
        if (dim(EXT.IN)[1]!=dim(weights)[1])
            stop("'EXT.IN' must have the same number of rows than 'weights'")
        if (!is.list(EXT.FUN))
            stop("'EXT.FUN' must be a list")
        if (!is.character(EXT.FUN[[1]]))
            stop("'EXT.FUN' must contain strings")
        if (dim(EXT.IN)[2]!=length(EXT.FUN))
            stop("'EXT.FUN' must contain a number of functions equal to the number of coloumns in 'EXT.IN'")
    }

    ris <- .Call("simulateprofiles",
        as.matrix(weights),
        as.list(Rules),
        as.character(f.pr.and),
        as.vector(x0),
        as.vector(Xmax),
        as.vector(lambda),
        as.character(act.fun),
        as.vector(alpha),
        as.vector(beta),
        as.vector(times),
        as.character(method),
        as.matrix(EXT.IN),
        as.list(EXT.FUN)
    )
    return(ris)
}

Dati <- matrix(c(1.793561,1.793561,1.793561,1.793561,1.793561,4.239327,4.239327,4.239327,4.239327,4.239327,1.031338,1.031338,1.031338,1.031338,1.031338,1.789676,1.789676,1.789676,1.789676,1.789676,7.767749,7.767749,7.767749,7.767749,7.767749), 5, 5)

pesi <- matrix(c(0,0,0,-1,0,-1,0,-1,0,0,0,-1,0,0,-1,-1,0,0,0,0,0,0,-1,0,0), 5, 5)

lambda <-c(1.1056752,0.9747519,1.0007163,1.0061211,1.1354204)

alpha <- c(9.869912,9.955962,10.056057,9.922354,9.950865)

beta <- c(0.4934634,0.4930697,0.4926102,0.4976155,0.4905647)

R <- NULL
f.pr.and <- NULL
N<-dim(Dati)[1]
X0<-runif(N,0,1)
Xmax<-rep(10,N)
ext.in <- NA
times<-seq(1,5)

ris <- simulateprofiles_dll(pesi, R, f.pr.and, X0, Xmax, lambda, "sigmoidal", alpha, beta, times, "Euler", ext.in, list())
ris
cat("ok")
