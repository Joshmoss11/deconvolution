library(nnls)
source('entropy_sites.R')
nnls_matrix <- function(X,Y){
	nans <- apply(X,1,function(x) {any(is.na(x))}) | apply(Y,1,function(x) {any(is.na(x))})
	X <- X[!nans,]
	Y <- as.matrix(Y[!nans,])
	
	N <- nrow(X); M <- ncol(X); K <- ncol(Y)

	nnls.betas <- matrix(NA,K,M)
	for (k in 1:K){
		nnls.betas[k,] <- nnls(X,Y[,k])$x
	}
	rownames(nnls.betas) <- colnames(Y)
	colnames(nnls.betas) <- colnames(X)
	return(nnls.betas)
}
norm_nnls_betas <- function(nnls.betas){
	nnls.betas.norm <- t(apply(nnls.betas,1,function(x) {x/sum(x)}))
	rownames(nnls.betas.norm) <- rownames(nnls.betas)
	colnames(nnls.betas.norm) <- colnames(nnls.betas)
	return(nnls.betas.norm)
}
run_nnls_matrix <- function(X.f,Y.f,out.f,out.norm.f=NULL,normbetas=T,rm.cols=c('Vasc_endothelium','Skin','Liver','Omentum'),nrows=-1,entropy_sites=F,num_entropy_sites=100){
	print(paste('Reading in reference table from',X.f))
	X <- as.matrix(read.csv(X.f,row.names=1,nrow=nrows,check.names=F))
	print('Read in reference table')
	X <- X[,!(colnames(X) %in% rm.cols)]
	Y <- as.matrix(read.csv(Y.f,row.names=1,nrow=nrows,check.names=F))
	print('Read in test samples')
	if (entropy_sites){
		print('Selecting best sites by entropy')
		sites <- entropy_sites(X,num_entropy_sites)
		X <- X[sites,]
		Y <- Y[sites,]
	}
	print('Performing nnls on all samples')
	nnls.betas <- nnls_matrix(X,Y)
	print('Done calculatign betas. Printing results')
	write.csv(nnls.betas,out.f,row.names=T,quote=F)
	if (normbetas){
		print('Calculating and printing normalized betas')
		write.csv(norm_nnls_betas(nnls.betas),out.norm.f,row.names=T,quote=F)
	}
	print('Done performings nnls')
}
