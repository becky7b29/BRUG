
rm(list=ls())

#### Read data ####

ver = "50"
type = "bD_2w"

numfiles = 100
p = 6
n = 100
nedge = p*(p-1)/2

setwd(paste("~/Paper/Gen_data/v4x_seed123/gendata_p",p,"_v",ver,"/Xk",sep=""))
Xk = read.table(paste("data_mygraph_n",n,"_1.txt",sep=""))
for( i in 2:numfiles ){
	Xk = rbind(Xk,read.table(paste("data_mygraph_n",n,"_",i,".txt",sep="")))
}
Xk = as.matrix(Xk)



#### Read covariates ####

setwd(paste("~/Paper/Gen_data/v4x_seed123/gendata_p",p,"_v",ver,sep=""))

X = read.table("xi_mygraph_n1.txt")
if( ver == 13 ){
	X[,4] = X[,3]==3
	X[,5] = X[,3]==4
	X[,6] = X[,3]==5
	X[,3] = X[,3]==2
}
X = cbind(rep(1,numfiles),as.matrix(X))

npar = dim(X)[2]


#### Read parameters and test PG Gibbs ####

setwd(paste("~/Paper/Gen_data/v4x_seed123/gendata_p",p,"_v",ver,"/adjacency",sep=""))
adj_True = list()
for( i in 1:numfiles ){
	adj_True[[i]] = as.matrix(read.table(paste("adj_mygraph_n1_",i,".txt",sep="")))
}

gij = matrix(rep(NA,nedge*numfiles),nrow=nedge)
for( l in 1:numfiles ){
	temp_adj = as.matrix(adj_True[[l]])
	gij[,l] = temp_adj[which(upper.tri(temp_adj))]
}


#### Analysis ####

library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)

setwd("~/Paper/Code/Rcpp")
#sourceCpp("Wrapper_bD_v3_parallel_v3_ud.cpp")
sourceCpp("Wrapper_bD_v3_parallel_v3_ud_nondiag_v1.cpp")
#sourceCpp("Wrapper_bD_v3_parallel_v3_ud_nondiag_v2.cpp")


Beta = matrix(rep(0,nedge*numfiles),nrow=nedge); 
C = diag(p); 


est_D = array(NA,c(p,p,numfiles))
for( l in 1:numfiles ){
	if( n == 1 ){
		est_D[,,l] = (Xk[l,])%o%(Xk[l,]);
	}else{
		est_D[,,l] = t(Xk[(n*(l-1)+1):(n*l),])%*%(Xk[(n*(l-1)+1):(n*l),]);
	}
}

b_prior = 5; D_prior = apply(est_D,1:2,mean)/n; 


burnin = 5000
if(type=="bD_2w"){ nmc=15000 }
if(type=="bD_6w"){ nmc=55000 }
if(type=="bD_10w"){ nmc=95000 }
if(type=="bD_20w"){ nmc=195000 }
init_burnin = 3000
init_nmc = 3000
inloop_burnin = 0
inloop_nmc = 1
PG_interval = 1
bD_interval = 1

(numcores = defaultNumThreads())
#setThreadOptions(numThreads = numcores-1)

t = proc.time()
test = BRUG( b_prior, D_prior, n, Xk, X, C, Beta, burnin, nmc, init_burnin, init_nmc, inloop_burnin, inloop_nmc, PG_interval, bD_interval)
#test = BRUG( b_prior, D_prior, n, Xk, X, C, Beta, burnin, nmc, init_burnin, init_nmc, inloop_burnin, inloop_nmc, PG_interval)
(t = proc.time()-t)


adj_temp = matrix(0,nrow=p,ncol=p)
adj_est_mymodel = list()

adj_vec = apply(test$adj[,,(burnin+1):(burnin+nmc+1)],c(1,2),mean)
gij_est_mymodel = (adj_vec[-(p:1)*(p:1+1)/2,]>0.5)+0


#### compute estimated FDR ####

PPI_bound = 0.5
PPI = adj_vec[-(p:1)*(p:1+1)/2,]

gij_est_mymodel = ((PPI>=(1-PPI_bound))+0)
edge_detect = sum((gij_est_mymodel==1)&(gij==1))
FP = sum((gij_est_mymodel==1)&(gij==0))

test$accept_b/(burnin+nmc)
test$accept_D/(burnin+nmc)

D_est = apply(test$D[,,(burnin+1):(burnin+nmc+1)],c(1,2),mean)
round(D_est,2)
(b_est = mean(test$b[(burnin+1):(burnin+nmc+1)]))
mean(diag(round(D_est,2)))

sum(gij_est_mymodel)
edge_detect
(estFDR = sum((1-PPI)*(((1-PPI)<=PPI_bound)+0))/sum(((1-PPI)<=PPI_bound)+0))
(trueFDR = FP/sum(gij_est_mymodel))
edge_detect/sum(gij)
sum(gij)
n

library(xtable)

output = as.data.frame(matrix(c(n,b_est,mean(diag(round(D_est,2))),estFDR*100,trueFDR*100,edge_detect/sum(gij)*100,sum((gij_est_mymodel==0)&(gij==0))/sum(gij==0)*100,sum(gij_est_mymodel==gij)/(nedge*numfiles)*100),nrow=1))
colnames(output) = c("Sample","b","D","eFDR","FDR","Sensitivity","Specificity",
"Overall Accuracy")

print(xtable(output), include.rownames=FALSE)



setwd(paste("~/Paper/cppSim_v",ver,"_",type,sep=""))

save.image(paste("cppSim_v",ver,"_n",n,".RData",sep=""))
