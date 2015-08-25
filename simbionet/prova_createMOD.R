setwd(".")

dyn.load("create_mod.dll")

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

ris <- createMOD(m=4, auto=F)
cat("ok")