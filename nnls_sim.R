nrows=-1
source('nnls_matrix_functions.R')
sim.id <- args<-commandArgs(trailingOnly=TRUE)[1]
deconv.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution'
matrices.dir <- file.path(deconv.dir,'matrices')
X.f <- file.path(matrices.dir,'ref_wi_cancer.csv')
sim.dir <- file.path(deconv.dir,'simulations')
Y.dir <- file.path(sim.dir,'y')
Y.f <- file.path(Y.dir,paste0('set',sim.id,'.csv'))
nnls.dir <- file.path(sim.dir,'nnls')
nnls.f <- file.path(nnls.dir,paste0('set',sim.id,'.csv'))
dir.create(nnls.dir,showWarnings=F)
print(paste('Reading in reference table from',X.f))
X <- as.matrix(read.csv(X.f,row.names=1,nrow=nrows,check.names=F))
rm.cols=c('Skin','Liver','Omentum')
X <- X[,!(colnames(X) %in% rm.cols)]
print('Read in reference table')
Y <- as.matrix(read.csv(Y.f,row.names=1,nrow=nrows,check.names=F))
X.2 <- merge(X,Y,by=0); rownames(X.2) <- X.2[,1]; X.2 <- X.2[,2:ncol(X.2)]
X <- as.matrix(X.2[,1:(ncol(X.2)-1)])
Y <- as.matrix(X.2[,ncol(X.2),drop=F])
nnls.betas <- nnls_matrix(X,Y)
print('Done calculating betas. Printing results')
write.csv(nnls.betas,nnls.f,row.names=F,quote=F)
