setwd(".")

dyn.load("hmm.so")
dyn.load("create_mod.so")

createMOD <- function(m=4, auto=FALSE)
{
    stopifnot(is.numeric(m))
    stopifnot(is.logical(auto))
    ris <- .Call("createMOD",
                 as.integer(m),
                 as.logical(auto)
					)
    return(ris)
}

HMM<-function(N=50, Cf.cl=0.3, gamma=2.2, DEGREE=NULL, INdegree=c("free", "out"), MODULES, prior.p.subnet=NULL, max.con=12, feedback=TRUE, zero.nodes.with.indegree0=TRUE, sepgraph=TRUE, r.tol=0.1, a.tol=1, iter=1)
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

    INdegree<-match.arg(INdegree)

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

    if (!is.logical(feedback))
        stop("'feedback' must be TRUE or FALSE")

    if ((zero.nodes.with.indegree0!=TRUE)&(zero.nodes.with.indegree0!=FALSE))
        stop("'zero.nodes.with.indegree0' must be TRUE or FALSE")

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

    ris <- .Call("HMM",
        as.integer(N),
        as.numeric(Cf.cl),
        as.numeric(gamma),
        as.vector(DEGREE),
        as.character(INdegree),
        as.list(MODULES),
        as.vector(prior.p.subnet),
        as.integer(max.con),
        as.logical(feedback),
        as.logical(zero.nodes.with.indegree0),
        as.logical(sepgraph),
        as.numeric(r.tol),
        as.numeric(a.tol),
        as.integer(iter)
    )

    return (ris)
}

MOD <- createMOD(m=2,auto=TRUE)
N <- 1000
L <- length(MOD)
p <- runif(L)
p <- p / sum(p)
deg <- c()
somma <- N
for (i in 1:N) {
    tmp <- runif(1, 0, (N + 1) / N)
    if (somma - tmp > 0) {
        somma <- somma - tmp
    }
    else tmp <- 0.0
    deg <- c(deg, tmp)
}
if (somma > 0) deg <- c(deg, somma) else deg <- c(deg, 0.0)
#cat("p = ", p)
#cat("deg = ", deg)
sum(deg)


fp <- file("stat_C.tsv", "w")
P <- 10
s <- 2 ^ seq(1, 9)

for (i in 1:P) {
    cat(file=fp, sprintf("\t%d", i))
}
cat(file=fp, "\tMedia\tDev. std.\n")

times <- c()
for (n in s) {
    cat(file=fp, sprintf("%d", n))
    for (i in 1:P) {
        t <- system.time(ris <- HMM(N=n,MODULES=MOD,Cf.cl=0, iter=1))
        cat(file=fp, sprintf("\t%3.3f", t[1]))
        times <- c(times, t[1])
    }
    cat(file=fp, sprintf("\t%3.3f\t%3.3f\n", mean(times), sd(times)))
}
close(fp)
