nrows=-1
deconv.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution'
matrices.dir <- file.path(deconv.dir,'matrices')
Y <- as.matrix(read.csv(file.path(matrices.dir,'test_samples.csv'),row.names=1,nrow=nrows,check.names=F))
plots.dir <- file.path(deconv.dir,'plots')
pdf(file.path(plots.dir,'density_plots.pdf'))
for(i in 1:ncol(Y)){
	plot(density(Y[,i]),main=colnames(Y)[i])
}
dev.off()
