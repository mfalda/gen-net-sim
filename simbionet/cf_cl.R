library(SimBioNeT)

MOD<-createMOD(m=4,auto=FALSE)
net2<-SBNT(N=20,MODULES=MOD,Cf.cl=0.1,sepgraph=FALSE)
net2