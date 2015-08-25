\name{SBNT}
\alias{SBNT}
\title{Gene network topology simulator}
\description{This program is a simulator that generates directed network topologies, so as to reproduce some important features of biological networks, according to the parameter settings given as input.}
\usage{
SBNT(N=50,Cf.cl=0.3,gamma=2.2,DEGREE=NULL,INdegree="free",MODULES,prior.p.subnet=NULL,max.con=12,feedback=TRUE,zero.nodes.with.indegree0=TRUE,sepgraph=TRUE,r.tol=0.1,a.tol=1,iter=1)

}

\arguments{
  \item{N}{Number of nodes in the network.}
  \item{Cf.cl}{Average clustering coefficient of each sub-network in the graph (see References for further details).Default=0.4, leading to an average clustering coefficient in the entire network ranging between 0.1 and 0.3 as observed by (Barabasi and Albert (1999)) in protein networks.}
  \item{gamma}{Parameter of the power law distribution of the degree of connectivity of the nodes in the graph.Default=2.2  equal to the average value observed by (Albert and Barabasi (2000)) in protein networks.}
  \item{DEGREE}{Probability distribution of the out-degree. Default=NULL: power-law distribution will be used. Otherwise a vector of length (N+1) should be provided with, in position k, the average number of nodes with desired degree (k-1)}
  \item{INdegree}{The in-degree distribution: either "free", i.e. it is not constrained to follow any distribution or "out", i.e. it is constrained to follow the same distribution of the out-degree (in this latter case the algorithm does not always converge)}
  \item{MODULES}{list of L modules obtained by running function "createMOD(m,auto=FALSE)" if self-regulation is prohibited; "createMOD(m,auto=TRUE)" if self-regulation is allowed"}
  \item{prior.p.subnet}{Vector of L elements containing the a priori probability of having a sub-network motif of type L. Default=equal probability}
  \item{max.con}{Maximum number of regulators that each node in the network can have. Default=12}
  \item{feedback}{If FALSE no modules with feedback are allowed. Default=TRUE}
  \item{zero.nodes.with.indegree0}{If TRUE (default) every node must have at least in.degree 1}
  \item{sepgraph}{If FALSE networks consisting of separate subgraph are prohibited. Default is TRUE.}
  \item{r.tol}{relative tolerance parameter for the power law distribution. Default=0.1 (See Details).}
  \item{a.tol}{absolute tolerance parameter for the power law distribution. Default=0.1 (See Details).}
  \item{iter}{number ofiterations. Default=1}
}

\details{}


\value{The function returns the adjacency matrix containing information on network connectivity: a value 1 in position (i,j) means that node j regulates
	     node i.
}

\author{
  Barbara Di Camillo
}

\references{
\item{1.}{Barbara Di Camillo, Marco Falda, Gianna Toffolo, Claudio Cobelli. Topology Simulation of Biological Networks:
a Hierarchical Modular Model. \emph{IEEE Transactions on Biomedical Engineering} \bold{submitted}.}
}


\seealso{\code{\link{UNSBNT}} generating undirected networks.}


\examples{
#MOD<-createMOD(m=4,auto=FALSE)
#net1<-SBNT(N=50,MODULES=MOD,sepgraph=FALSE)
#net2<-SBNT(N=50,MODULES=MOD,Cf.cl=0,sepgraph=FALSE)
}



\keyword{datagen}

