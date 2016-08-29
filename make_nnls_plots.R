library(nnls)
library(gplots)
deconv.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution'
nnls.dir <- file.path(deconv.dir,'nnls')
plots.dir <- file.path(deconv.dir,'plots')
dir.create(plots.dir)
x <- read.csv(file.path(nnls.dir,'nnls_betas_norm.csv'),row.names=1,check.names=F)
pdf(file.path(plots.dir,'betas_heatmap.pdf'))
heatmap.2(as.matrix(x),trace='none',margins = c(8, 8))
dev.off()

pal <- rainbow(length(lab))
pdf(file.path(plots.dir,'betas_pies.pdf'))
x[x<0.005] <- 0
lab <- c(colnames(x),'Other')
for (i in 1:nrow(x)){
	x.cur <- x[i,]
	pie(x.cur[x.cur!=0],labels=lab[x.cur!=0],col=pal[x.cur!=0],main=rownames(x)[i])
}
dev.off()
