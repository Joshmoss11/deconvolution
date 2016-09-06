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
