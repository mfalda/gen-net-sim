\name{simulatenet}
\alias{simulatenet}
\title{Gene network simulator integrating Boolean logic in a continuous dynamic model of regulation}
\description{This program is a simulator that generates network topology, regulatory rules and expression profiles so as to reproduce some important features of regulation in real gene networks, according to the parameter settings given as input.}
\usage{
simulatenet <- function(N=50, connectivity="MTM", max.reg=12, gamma=2.2, INdegree="free", Cf.cl=0.4, num.subnet=c(5,5,10), kappa=3, f.pr.and="", act.fun="sigmoidal", alpha=c(10,0.2), beta=c(0.5,0.01), lambda=c(1,0.1), Xmin=NULL, Xmax=NULL, X0=c(0,1), weight.par=c(1,0), params=c(2,2,2,0,0,1), times=seq(0,5,0.05), stat_thr=0.001, stat_width=0, method="rkf45", itera=1, save=FALSE)
}

\arguments{
  \item{N}{Number of genes in the network.}
  \item{connectivity}{Type of connectivity in the network. One among "random","scale free", "MTM" and "geometric". Default="MTM"}
  \item{max.reg}{Maximum number of regulators that each node (gene) in the network can have. Default=12}
  \item{gamma}{In case "scale-free" or "MTM" topology is chosen, this is the parameter of the power law distribution of the degree of connectivity of the nodes in the graph.Default=2.2  equal to the average value observed by (Albert and Barabasi (2000)) in protein networks.}
  \item{INdegree}{The in-degree distribution (applicable only if connectivity="MTM"): either "free", i.e. it is not constrained to follow any distribution or "out", i.e. it is constrained to follow the same distribution of the out-degree (in this latter case the algorithm does not always converge)}
  \item{Cf.cl}{Average clustering coefficient of each sub-network in the graph (see References for further details).Default=0.4, leading to an average clustering coefficient in the entire network ranging between 0.1 and 0.3 as observed by (Barabasi and Albert (1999)) in protein networks.}
  \item{num.subnet}{Vector of 3 elements containing the maximum number of nodes in each sub-network motif of type 1, 2 or 3 respectively (see References for further details). Default=c(5,5,10)}
  \item{kappa}{average number of regulators that each node (gene) in the network has with random topology. Default=3}
  \item{f.pr.and}{Function with domain and codomain in [0-1] that expresses the probability to obtain a cooperative rather than a synergic rule as a function of the level of the binary tree used to implement the regulatory interactions among the regulators. Given a binary tree of maximum L levels, the probability to have a cooperative rule at level i is f.pr.and(i/L), whereas the probability to have a synergic rule is 1-f.pr.and(i/L). See Reference for further details}
  \item{act.fun}{The activation function: either "linear", i.e. the target function resulting from regulatory rules or "sigmoidal", i.e. the modulated target function (see References for further details). Default="sigmoidal"}
  \item{alpha}{Vector of parameters of the Activation sigmoid function. If NULL, alpha is set equal to rnorm(1,10,0.2) for each gene.}
  \item{beta}{Vector of parameters of the Activation sigmoid function. If NULL, alpha is set equal to rnorm(1,0.5,0.01) for each gene.}
  \item{lambda}{vector of time constants influencing both the rate of transcription and the spontaneous degradation term for each gene i. If NULL, lambda is set equal to rnorm(N,1,0.1).}
  \item{Xmin}{vector of minimum level of expression of genes. If NULL, Xmin is set equal to rep(0,N).}
  \item{Xmax}{vector of maximum level of expression of genes. If NULL, Xmax is set equal to rep(10,N).}
  \item{X0}{Initial conditions, i.e. gene expression values at time 0 scaled between 0 and 1. If NULL, X0 is set equal to runif(1,0,1) for each node i.}
  \item{weight.par}{vector of two elements: mean and sd of the Gaussian distribution used to sample regulatory efficiency w$_{ir}$ of regulator r on gene i (see References for further details). Default=c(1,0.2)}
  \item{params}{A vector identifying the type of distribution for input parameters lambda, alpha, beta, Xmin, Xmax and X0; 0 stands for "external", 1 for "uniform", 2 for "normal" and 3 for "log-normal"; the last element can assume only values 0 or 1. Default=c(2,2,2,0,0,1)}
  \item{times}{Time samples at which explicit estimates of gene expression are desired. The first value must be 0. Default=seq(0,5,0.05)}
   \item{stat_thr}{the absolute difference between consecutive points in the generated profiles that, when satisfied by all the profiles, allows to consider the solution in a stationary point. Default=0.001}
   \item{stat_width}{the percentage of consecutive sampling times during which the "stat_thr" threshold has to be satisfied in order to consider the solution in a stationary point. Default=0}
  \item{method}{One among "Euler", "rkf45", "rkrk", "rk8pd", "rk2imp", "rk4imp", "gear1", "gear2". Method used to solve differential equations. Either the Euler method or one of the GNU Scientific Library ODE solvers. Default="rkf45"}
  \item{itera}{number of iterations. Default=1}
  \item{save}{If TRUE the connectivity matrix, the generated expression profiles, the regulatory rules and a matrix with parameters lambda, alpha, beta, Xmin, Xmax and X0 are saved. Default=FALSE}
}

\details{

The algorithm, implemented in R, requires the user to specify the number N of genes
in the network. For all other parameters defaults are available, but the user can set different
values, if desired.
The simulator first generates the topology, to assign regulators to each gene, with a constraint
on the maximum number of regulators that each gene can have. With random topology the parameter kappa is used to set the average in-degree of each node; with scale-free topology, the parameter gamma of the power-law distribution has to be set; with MTM topology, besides gamma, average clustering coefficient Cl.cf and maximum number of nodes in each sub-network structure are required.
Once the topology is given, regulators of each node i are assigned. Based on them, a regulatory
rule is simulated by a random combination of AND OR NOT rules, extended to fuzzy logic. efficiency parameters w$_{ir}$ (r=1, ..., Ri) are randomly sampled from a Gaussian distribution to account for different regulation efficiency of the regulators of each gene i (absolute values are taken in case of negative values).
Once the network is completely specified both in terms of topology and rules, expression profiles
are generated using differential equations (see References for further details).
Data are generated as continuous variables in the range Xmin$_{i}$-Xmax$_{i}$ (for each gene i) observed with
a sampling grid given by the user.

The possibility to generate gene expression data given topology, regulatory rules and kinetic
parameters is implemented in function "simulateprofiles".


}


\value{The function returns two lists, the first element is a list with five elements:
  \enumerate{
  \item{expr.data: the matrix of gene expression data with genes in rows and time-samples in columns;}
  \item{weight.matrix: the adjacence matrix containing information on network connectivity and regulatory
	     efficiency and sign: a non zero value in position (i,j) means that gene j regulates
	     gene i with efficiency and sign indicated by the value itself;}
   \item{lambda: a vector of N elements containing parameters lambda setted for each gene.}
   \item{alpha: a vector of N elements containing parameters alpha setted for each gene.}
   \item{beta: a vector of N elements containing parameters beta setted for each gene.}
   }
  The second element contains a list with the regulatory Rules. The regulatory rule for gene i is codified in position [[i]] of the list as a vector representing a structure of a binary tree (see help(simulateprofiles) for further details.
}


\author{
  Barbara Di Camillo and Marco Falda
}

\references{
\enumerate{
\item{Barbara Di Camillo, Gianna Toffolo, Claudio Cobelli. A Gene Network Simulator to Assess Reverse Engineering Algorithms.
  \emph{NYAS} \bold{Submitted}.}
\item{Barbara Di Camillo, Gianna Toffolo, Claudio Cobelli. A gene network simulator integrating Boolean regulation in a continuous dynamic model
  \emph{5th European Conference on Computational Biology - ECCB '06} \bold{Eilat, Israel, January 20-24, 2007}.}
  }
}


\seealso{\code{\link{simulateprofiles}} implementing the possibility to generate gene expression
data given topology, regulatory rules and kinetic parameters.}


\examples{
net<-simulatenet(N=50,connectivity="MTM",act.fun="sigmoidal", method="rkf45")
Data<-net[[1]][[1]]
W.matrix<-net[[1]][[2]]
R<-net[[2]]
param.lambda<-net[[1]][[3]]
param.alpha<-net[[1]][[4]]
param.beta<-net[[1]][[5]]

net<-simulatenet(N=10,connectivity="random",kappa=5, method="rkf45")

net<-simulatenet(N=50,connectivity="MTM",gamma=2.5, Cf.cl=0.2, act.fun="sigmoidal", method="rkf45")

net <- simulatenet(N=50, ,connectivity="MTM", f.pr.and="1 / (1 + _e ^ (-10 * (x - 0.5)))", act.fun="sigmoidal", method="rkf45")

}



\keyword{datagen}
\keyword{ts}

