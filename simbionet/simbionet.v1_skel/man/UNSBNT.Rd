\name{UNSBNT}
\alias{UNSBNT}
\title{Gene network topology simulator}
\description{This program is a simulator that generates directed network topologies, so as to reproduce some important features of biological networks, according to the parameter settings given as input.}
\usage{
UNSBNT(N=50,Cf.cl=0.3,gamma=2.2,DEGREE=NULL,INdegree="free",MODULES,prior.p.subnet=NULL,max.con=12,feedback=TRUE,zero.nodes.with.indegree0=TRUE,sepgraph=TRUE,r.tol=0.1,a.tol=1,iter=1)

}

\arguments{
  \item{N}{Number of nodes in the network.}
  \item{Cf.cl}{Average clustering coefficient of each sub-network in the graph (see References for further details).Default=0.4, leading to an average clustering coefficient in the entire network ranging between 0.1 and 0.3 as observed by (Barabasi and Albert (1999)) in protein networks.}
  \item{gamma}{Parameter of the power law distribution of the degree of connectivity of the nodes in the graph.Default=2.2  equal to the average value observed by (Albert and Barabasi (2000)) in protein networks.}
  \item{DEGREE}{Probability distribution of the degree. Default=NULL: power-law distribution will be used. Otherwise a vector of length (N+1) should be provided with, in position k, the average number of nodes with desired degree (k-1)}
  \item{MODULES}{list of L modules obtained by running function "createMOD(m=4,auto=FALSE)" if self-regulation is prohibited; "createMOD(m=4,auto=TRUE)" if self-regulation is allowed"}
  \item{prior.p.subnet}{Vector of L elements containing the a priori probability of having a sub-network motif of type L. Default=equal probability}
  \item{max.con}{Maximum degree that each node in the network can have. Default=12}
  \item{sepgraph}{If FALSE networks consisting of separate subgraph are prohibited. Default is TRUE.}
  \item{r.tol}{relative tolerance parameter for the power law distribution. Default=0.1 (See Details).}
  \item{a.tol}{absolute tolerance parameter for the power law distribution. Default=0.1 (See Details).}
  \item{iter}{number of iterations. Default=1}
}

\details{}


\value{The function returns the adjacency matrix containing information on network connectivity. The matrix is symmetric.
}

\author{
  Barbara Di Camillo
}

\references{
\item{1.}{Barbara Di Camillo, Marco Falda, Gianna Toffolo, Claudio Cobelli. Topology Simulation of Biological Networks:
a Hierarchical Modular Model. \emph{IEEE Transactions on Biomedical Engineering} \bold{submitted}.}
}


\seealso{\code{\link{SBNT}} generating directed networks.}


\examples{
#MOD<-createMOD(m=4,auto=FALSE)
#net1<-UNSBNT(N=50,MODULES=MOD,sepgraph=FALSE)
#net2<-UNSBNT(N=50,MODULES=MOD,Cf.cl=0,sepgraph=FALSE)
}



\keyword{datagen}

