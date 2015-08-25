simulateprofiles<-function(weights, Rules=list(), f.pr.and="", x0=c(0,1), Xmin=NULL, Xmax=NULL, lambda=c(1,0.1), act.fun="sigmoidal", alpha=c(10,0.2), beta=c(0.5,0.01), times=seq(0,5,0.05), method="rkf45", EXT.IN=NA, EXT.FUN=list(), params=c(2,2,2,0,0,1), ko.experim=NULL, save=FALSE, itera=1)
{
    if (!is.null(weights)) {
        if (!is.numeric(weights))
            stop("'weights' must be numeric")
        if (!is.matrix(weights))
            stop("'weights' must be a matrix")
        N<-dim(weights)[1]
    }
    else
        N<-dim(read.table("weights1.txt"))

    if (!is.null(Rules)) {
        if (!is.list(Rules))
            stop("'rules' must be a list or NULL")
        if (length(Rules) > 0 && !is.numeric(Rules[[1]]))
            stop("'Rules' must contain numeric values")
    }

    if (!is.null(f.pr.and)) {
        if (!is.character(f.pr.and))
            stop("`f.pr.and' must be NULL or a string")
    }
    else
        f.pr.and <- ""

    if (!is.numeric(lambda))
        stop("`lambda' must be numeric")

    if (!is.numeric(x0))
        stop("'x0' must be numeric")

    if (is.null(Xmin))
        Xmin<-rep(0,N)
    else if (!is.numeric(Xmin))
        stop("`Xmin' must be numeric")

    if (is.null(Xmax))
        Xmax<-rep(10,N)
    if (!is.numeric(Xmax))
        stop("'Xmax' must be numeric")

    if (!is.numeric(lambda))
        stop("'lambda' must be numeric")

    if (!is.numeric(alpha))
        stop("`alpha' must be numeric")

    if (!is.numeric(beta))
        stop("`beta' must be numeric")

    if (!is.numeric(times))
        stop("'times' must be numeric")

    if (!is.list(EXT.FUN))
        stop("'EXT.FUN' must be a list")

    if (length(EXT.FUN)>0) {
        if (!is.numeric(EXT.IN))
            stop("'EXT.IN' must be numeric")
        if (!is.matrix(EXT.IN))
            stop("'EXT.IN' must be a matrix")
        if (!is.character(EXT.FUN[[1]]))
            stop("'EXT.FUN' must contain strings")
        if (dim(EXT.IN)[2]!=length(EXT.FUN))
            stop("'EXT.FUN' must contain a number of functions equal to the number of coloumns in 'EXT.IN'")
    }

    if (!is.numeric(params))
        stop("`params' must be numeric")

    if (!is.null(ko.experim) && !is.numeric(ko.experim))
        stop("`ko.experim' must be NULL or numeric")

    if (!is.logical(save))
        stop("`save' must be logical")

    if (!is.numeric(itera))
        stop("`itera' must be numeric")

    ris <- .Call("simulateprofiles",
        as.integer(N),
        weights,
        Rules,
        as.character(f.pr.and),
        as.vector(x0),
        as.vector(Xmin),
        as.vector(Xmax),
        as.vector(lambda),
        as.character(act.fun),
        as.vector(alpha),
        as.vector(beta),
        as.vector(times),
        as.character(method),
        as.matrix(EXT.IN),
        as.list(EXT.FUN),
        as.vector(params),
        as.vector(ko.experim),
        as.vector(save),
        as.integer(itera)
    )
    return(ris)
}
