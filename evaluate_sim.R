sim.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution/simulations'
betas.actual.dir <- file.path(sim.dir,'betas')
groups.actual.dir <- file.path(sim.dir,'groups')
betas.nnls.dir <- file.path(sim.dir,'nnls')
betas.matrices.dir <- file.path(sim.dir,'betas_matrices')
dir.create(betas.matrices.dir,showWarnings=F)
n.sim <- length(list.files(betas.nnls.dir))
betas.nnls.1 <- read.csv(file.path(betas.nnls.dir,paste0('set','1','.csv')),check.names=F)
n.tissues <- ncol(betas.nnls.1)
tissues <- colnames(betas.nnls.1)
betas.nnls <- matrix(0,n.sim,n.tissues)
betas.actual <- matrix(0,n.sim,n.tissues)

for (i in 1:n.sim){
	betas.nnls.f <- file.path(betas.nnls.dir,paste0('set',i,'.csv'))
	betas.nnls[i,] <- as.numeric(read.csv(betas.nnls.f))
	betas.actual.f <- file.path(betas.actual.dir,paste0('set',i,'.csv'))
	groups.actual.f <- file.path(groups.actual.dir,paste0('set',i,'.csv'))
	groups.actual <- read.csv(groups.actual.f,header=F,stringsAsFactors=F)[,1]
	groups.idx <- match(groups.actual,tissues)
	betas.actual.i <- read.csv(betas.actual.f,header=F,stringsAsFactors=F)[,1]
	betas.actual[i,groups.idx] <- betas.actual.i
}

colnames(betas.actual) <- tissues
colnames(betas.nnls) <- tissues
write.csv(betas.actual,file.path(betas.matrices.dir,'betas_actual.csv'),row.names=F)
write.csv(betas.nnls,file.path(betas.matrices.dir,'betas_nnls.csv'),row.names=F)

cors <- sapply(1:nrow(betas.nnls),function(i) cor(betas.actual[i,],betas.nnls[i,]))
pdf(file.path(sim.dir,'simulation_eval.pdf'))
for (i in 1:nrow(betas.actual)){
    plot(betas.actual[i,],betas.nnls[i,],xlab='Simulated betas',ylab='Predicted betas',main=paste0('Set',i))
    text(betas.actual[i,],betas.nnls[i,],labels=colnames(betas.actual))
}
dev.off()
