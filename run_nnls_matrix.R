nrows=-1
library(nnls)
nnls_matrix <- function(X,Y){
	nans <- apply(X,1,function(x) {any(is.na(x))}) | apply(Y,1,function(x) {any(is.na(x))})
	X <- X[!nans,]
	Y <- Y[!nans,]
	
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
	rownames(nnls.betas.norm) <- colnames(Y)
	colnames(nnls.betas.norm) <- colnames(X)
}
run_nnls_matrix <- function(X.f,Y.f,out.f,out.norm.f=NULL,normbetas=T,rm.cols=c('Vasc_endothelium','Skin',
												'Liver','Omentum')){
	X <- as.matrix(read.csv(X.f,row.names=1,nrow=nrows,check.names=F))[,!(colnames(X) %in% rm.cols)]
	Y <- as.matrix(read.csv(Y.f,row.names=1,nrow=nrows,check.names=F))
	nnls.betas <- nnls_matrix(X,Y)
	write.csv(nnls.betas,out.f,row.names=T,quote=F)
	if (normbetas){
		write.csv(norm_nnls_betas(nnls.betas),out.norm.f,row.names=T,quote=F)
	}
}

deconv.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution'
matrices.dir <- file.path(deconv.dir,'matrices')
nnls.dir <- file.path(deconv.dir,'nnls')
dir.create(nnls.dir,showWarnings=F)
X.f <- file.path(matrices.dir,'ref_no_cancer.csv')
Y.f <- file.path(matrices.dir,'test_samples.csv')
out.f <- file.path(nnls.dir,'nnls_betas.csv')
out.norm.f <- file.path(nnls.dir,'nnls_betas_norm.csv')
X.f.w.c <- file.path(matrices.dir,'ref_wi_cancer.csv')
out.f.w.c <- file.path(nnls.dir,'nnls_betas_wi_cancer.csv')
out.norm.f.w.c <- file.path(nnls.dir,'nnls_betas_norm_wi_cancer.csv')

run_nnls_matrix(X.f,Y.f,out.f,out.norm.f)
run_nnls_matrix(X.f.wc,Y.f,out.f.wc,out.norm.f.wc)
