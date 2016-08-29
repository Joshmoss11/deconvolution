nrows=-1
deconv.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution'
matrices.dir <- file.path(deconv.dir,'matrices')
X <- as.matrix(read.csv(file.path(matrices.dir,'ref_no_cancer.csv'),row.names=1,nrow=nrows,check.names=F))
d <- dist(t(X))
pdf(file.path(deconv.dir,'plots','ref_tissues_cluster.pdf'))
plot(hclust(d))
dev.off()
