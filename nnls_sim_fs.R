library(FSelector)
nrows=1000
source('nnls_matrix_functions.R')
rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
get_max_diff <- function(mat){
    mat.min <- apply(mat,1,min)
    mat.max <- apply(mat,1,max)
    return(mat.max-mat.min)
}

nnls_matrix2 <- function(X,Y){
    nnls.betas <- apply(Y,2,function(y) nnls(X,y)$x)
    return(t(nnls.betas))
}
deconv.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution'
matrices.dir <- file.path(deconv.dir,'matrices')
X.f <- file.path(matrices.dir,'ref_wi_cancer.csv')
sim.dir <- file.path(deconv.dir,'simulations')
fs.dir <- file.path(sim.dir,'fs')
dir.create(fs.dir,showWarnings=F)
Y.dir <- file.path(sim.dir,'y')
betas.actual.dir <- file.path(sim.dir,'betas')
groups.actual.dir <- file.path(sim.dir,'groups')
print(paste('Reading in reference table from',X.f))
X <- as.matrix(read.csv(X.f,row.names=1,nrow=nrows,check.names=F))
rm.cols=c('Skin','Liver','Omentum')
X <- X[,!(colnames(X) %in% rm.cols)]
print('Read in reference table')
n.sim <- length(list.files(Y.dir))
n.tissues <- ncol(X)
tissues <- colnames(X)
betas.actual <- matrix(0,n.sim,n.tissues)
for (i in 1:n.sim){
	Y.f <- file.path(Y.dir,paste0('set',i,'.csv'))
	Y <- as.matrix(read.csv(Y.f,row.names=1,nrow=nrows,check.names=F,header=F))
	X.2 <- merge(X,Y,by=0); rownames(X.2) <- X.2[,1]; X.2 <- X.2[,2:ncol(X.2)]
	X <- as.matrix(X.2[,1:(ncol(X.2)-1)])
	Y <- as.matrix(X.2[,ncol(X.2),drop=F])
	colnames(Y)[1] <- paste0('set',i)
	if (i==1){ Y.all <- Y
	} else { Y.all <- merge(Y.all,Y,by=0); rownames(Y.all) <- Y.all[,1]; Y.all <- Y.all[,2:ncol(Y.all)]
	}
    betas.actual.f <- file.path(betas.actual.dir,paste0('set',i,'.csv'))
    groups.actual.f <- file.path(groups.actual.dir,paste0('set',i,'.csv'))
    groups.actual <- read.csv(groups.actual.f,header=F,stringsAsFactors=F)[,1]
    groups.idx <- match(groups.actual,tissues)
    betas.actual.i <- read.csv(betas.actual.f,header=F,stringsAsFactors=F)[,1]
    betas.actual[i,groups.idx] <- betas.actual.i
}

nans <- apply(X,1,function(x) {any(is.na(x))}) | apply(Y.all,1,function(x) {any(is.na(x))})
X <- X[!nans,]
Y.all <- as.matrix(Y.all[!nans,])

evaluator <- function(subset){
	if (length(subset)==1 && subset==start){
		return(0.01)
	} else if (length(subset)==1 && subset!=start){
		return(0)
	} else {
		betas.nnls <- nnls_matrix2(X[subset,],Y.all[subset,])
		cors <- sapply(1:nrow(betas.nnls),function(i) cor(betas.actual[i,],betas.nnls[i,]))
  		print(subset)
  		print(mean(cors))	
		return(mean(cors))
	}
}
X.max.diff <- get_max_diff(X)
keep <- X.max.diff>0.3
X <- X[keep,]
Y.all <- Y.all[keep,]

cgs <- rownames(X)
start <- sample(cgs,1)
cgs.subset <- forward.search(cgs,evaluator)
write.table(cgs.subset,file.path(fs.dir,'cgs.csv'),sep=',',row.names=F,col.names=F,quote=F)
#print('Done calculating betas. Printing results')
#write.csv(nnls.betas,nnls.f,row.names=F,quote=F)
