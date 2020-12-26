
##### keep sparsity and change only 1 edge for adjacent graphs #####

rm(list=ls())
set.seed(123)

## initialize paramaters

numfiles <- 100			# number of groups
p <- 6					# number of variables estimated
nedge <- p*(p-1)/2

pij <- rep(NA,nedge)	# vector for storing p(gij=1)'s
cij <- rep(NA,nedge)	# vector for storing upper triangular entries of adjacency matrix
aij <- rep(NA,nedge)	# vector for storing alpha_ij's

for(k in 1:numfiles) { 
    nam <- paste("adj", k, sep = "")
    assign(nam, matrix(rep(0,p*p),nrow=p))
}

npar <- 1
x <- matrix(seq(0,50,length.out=numfiles),ncol=npar)	# v39
#x <- matrix(seq(0,50*2,length.out=numfiles),ncol=npar)	# v39_1 wider x range
#x <- matrix(seq(0,125,length.out=numfiles),ncol=npar)	# v39_250_1
bij <- matrix(rep(NA,nedge*npar),ncol = npar)		# vector for storing beta_ij's


#### create adjacency matrices ####

aij <- matrix(c(-60,rep(-2.1,14)),ncol=1) # overall sparsity
bij[,1] <- matrix(c(2,rep(0,14)),ncol=1) # across-group sparsity

pij <- matrix(rep(NA,nedge*numfiles),ncol = numfiles)
gij <- matrix(rep(NA,nedge*numfiles),ncol = numfiles)

for(i in 1:numfiles){
	
	adj <- matrix(rep(0,p*p),nrow=p)
	
	pij[,i] <- exp(aij+apply(t(bij)*x[i,],2,sum))/(1+exp(aij+apply(t(bij)*x[i,],2,sum)))
	cij <- sapply((1:nedge),function(x){sample(c(0,1),size=1,prob=c(1-pij[x,i],pij[x,i]))})
	adj[which(upper.tri(adj))] <- cij
	
	assign(paste("adj",i,sep=""),adj)
	for( j in 1:nedge){ gij[j,i] <- adj[which(upper.tri(adj))][j] }
	
}




for(i in 1:numfiles){
	
	adj <- matrix(rep(0,p*p),nrow=p)
	
	adj[which(upper.tri(adj))] <- c(gij[1,i],gij[2:nedge,1])
	adj = adj + t(adj) + diag(p)
	
	assign(paste("adj",i,sep=""),adj)
	
}

gij[2:nedge,2:numfiles] = gij[2:nedge,1]

apply(gij,1,sum)
apply(gij,2,sum)
sum(gij)



## function for generating multivariate normal data

gendata <- function(adj,n,graphno){

	nk <- n		# sample size of group k
	p <- dim(adj)[1]
	
	Xk <- matrix(rep(NA,nk*p),nrow=nk)	# matrix for storing the sample group k
	Zk <- matrix(rep(NA,p*nk),ncol=p)	# matrix for storing N(0,1) r.v.
	
	C <- matrix(as.matrix(rgwish(1,adj,b=p+10,D=diag(p))),ncol=p) 
	Sigma <- solve(C)
	
	if(sum(eigen(C)$values>0)==p){	# ensure the precision matrix is positive definite
		
		# eigendecomposition of Sigma where columns of U are unit vectors
#		D <- diag(eigen(Sigma)$values)
#		U <- sapply(1:p,function(x){eigen(Sigma)$vectors[,x]/norm(as.matrix(eigen(Sigma)$vectors[,x]^2))})
		
#		Zk <- matrix(rnorm(p*nk),ncol=p)
#		Xk <- t(U%*%sqrt(D)%*%t(Zk))
		Xk = matrix(mvrnorm(n=nk,rep(0,p),Sigma),nrow=nk)
	}else{
		warning("The precision matrix is not positive definite")
	}
	
	setwd(paste("F:\\Paper\\Gen_data\\v4x_seed123\\gendata_p",p,"_v41",sep=""))
	write.matrix(Xk,file=paste("Xk\\data_mygraph_n",n,"_",graphno,".txt",sep=""),sep=" ",nk)
	write.matrix(adj,file=paste("adjacency\\adj_mygraph_n",n,"_",graphno,".txt",sep=""),sep=" ",p)
	write.matrix(C,file=paste("Omega\\Omega_mygraph_n",n,"_",graphno,".txt",sep=""),sep=" ",p)
	write.matrix(Sigma,file=paste("Sigma\\Sigma_mygraph_n",n,"_",graphno,".txt",sep=""),sep=" ",p)
	
}


## generate data

library(BDgraph)
library(MASS)

size = 1

for(j in 1:numfiles){
	gendata(get(paste("adj",j,sep="")),size,j)
}


setwd(paste("F:\\Paper\\Gen_data\\v4x_seed123\\gendata_p",p,"_v41",sep=""))
write.matrix(x,file=paste("xi_mygraph_n",size,".txt",sep=""),sep=" ",numfiles)
write.matrix(aij,file=paste("aij_mygraph_n",size,".txt",sep=""),sep=" ",numfiles)
write.matrix(bij,file=paste("bij_mygraph_n",size,".txt",sep=""),sep=" ",numfiles)










































