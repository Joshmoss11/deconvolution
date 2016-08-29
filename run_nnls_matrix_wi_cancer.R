nrows=-1
library(nnls)
deconv.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution'
matrices.dir <- file.path(deconv.dir,'matrices')
nnls.dir <- file.path(deconv.dir,'nnls')
dir.create(nnls.dir)
X <- as.matrix(read.csv(file.path(matrices.dir,'ref_wi_cancer.csv'),row.names=1,nrow=nrows,check.names=F))
Y <- as.matrix(read.csv(file.path(matrices.dir,'test_samples.csv'),row.names=1,nrow=nrows,check.names=F))

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
write.csv(nnls.betas,file.path(nnls.dir,'nnls_betas_wi_cancer.csv'),row.names=T,quote=F)

nnls.betas.norm <- t(apply(nnls.betas,1,function(x) {x/sum(x)}))
rownames(nnls.betas.norm) <- colnames(Y)
colnames(nnls.betas.norm) <- colnames(X)
write.csv(nnls.betas.norm,file.path(nnls.dir,'nnls_betas_norm_wi_cancer.csv'),row.names=T,quote=F)
