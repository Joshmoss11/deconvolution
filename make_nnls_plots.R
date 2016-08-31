library(nnls)
library(gplots)
make_nnls_pies <- function(matrix.f,out.f){
	x <- read.csv(matrix.f,row.names=1,check.names=F)
	lab <- c(colnames(x),'Other')
	pal <- rainbow(length(lab))
	pdf(out.f)
	x[x<0.005] <- 0
	for (i in 1:nrow(x)){
    	x.cur <- x[i,]
    	pie(x.cur[x.cur!=0],labels=lab[x.cur!=0],col=pal[x.cur!=0],main=rownames(x)[i])
	}
	dev.off()
}
make_nnls_heatmap <- function(matrix.f,out.f){
    x <- read.csv(matrix.f,row.names=1,check.names=F)
	pdf(out.f)
	heatmap.2(as.matrix(x),trace='none',margins = c(8, 8))
	dev.off()
}

deconv.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution'
nnls.dir <- file.path(deconv.dir,'nnls')
plots.dir <- file.path(deconv.dir,'plots')
dir.create(plots.dir,showWarnings=F)
matrix.f <- file.path(nnls.dir,'nnls_betas_norm.csv')
pies.f <- file.path(plots.dir,'betas_pies.pdf')
heatmap.f <- file.path(plots.dir,'betas_heatmap.pdf')

matrix.f.w.c <- file.path(nnls.dir,'nnls_betas_norm_wi_cancer.csv')
pies.f.w.c <- file.path(plots.dir,'betas_pies_wi_cancer.pdf')
heatmap.f.w.c <- file.path(plots.dir,'betas_heatmap_wi_cancer.pdf')

make_nnls_pies(matrix.f,pies.f)
make_nnls_heatmap(matrix.f,heatmap.f)

make_nnls_pies(matrix.f.w.c,pies.f.w.c)
make_nnls_heatmap(matrix.f.w.c,heatmap.f.w.c)
