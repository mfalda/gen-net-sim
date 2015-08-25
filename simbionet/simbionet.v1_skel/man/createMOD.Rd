\name{createMOD}
\alias{createMOD}
\title{Gene network topology simulator}
\description{This program is a simulator that generates directed network topologies, so as to reproduce some important features of biological networks, according to the parameter settings given as input.}
\usage{
createMOD(m=4,auto=FALSE)
}

\arguments{
  \item{m}{maximum size of constitutive modules. Default m=4}
  \item{auto}{if auto=FALSE modules with autoregulation are not allowed. Default auto=FALSE}
  }

\details{}


\value{The function upload the set of modules useful to call SBNT and UNSBNT.
}

\author{
  Barbara Di Camillo
}

\references{
\item{1.}{Barbara Di Camillo, Marco Falda, Gianna Toffolo, Claudio Cobelli. Topology Simulation of Biological Networks: 
a Hierarchical Modular Model. \emph{IEEE Transactions on Biomedical Engineering} \bold{submitted}.}
}


\seealso{}


\examples{
MOD<-createMOD(m=4,auto=FALSE)
}



\keyword{datagen}

