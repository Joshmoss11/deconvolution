library(nnls)
rows=-1

data.450k <- read.csv('betas_means.csv',header=T,row.names=1,nrow=rows,check.names=F)
data.epic <- read.csv('betas_epic.csv',header=T,row.names=1,nrow=rows,check.names=F)
data.epic.2 <- read.csv('processed_betas/MI_tech.csv',header=T,row.names=1,nrow=rows,check.names=F)
data.epic <- merge(data.epic,data.epic.2,by=0)
rownames(data.epic) <- data.epic[,1]; data.epic <- data.epic[,2:ncol(data.epic)]


data.all <- merge(data.450k,data.epic,by=0)
rownames(data.all) <- data.all[,1]; data.all <- data.all[,2:ncol(data.all)]
p.idx <- c(25:28,73:88)
#ref.idx <- c(4:12,17:23,47,50:67,71:72)
ref.idx <- c(4:12,14,17:21,23,29:45,47,51:67,71:72)
#Y <- as.matrix(data.all[,p.idx])
X <- as.matrix(data.all[,ref.idx])

Y.sim <- read.csv('simulation_y.csv',row.names=1)
m <- match(rownames(X),rownames(Y.sim))
Y.sim <- Y.sim[m,]
nnls.sim <- matrix(nrow=ncol(Y.sim),ncol=ncol(X))
for (i in 1:ncol(Y.sim)){
	nnls.sim[i,] <- nnls(X,Y.sim[,i])$x
	nnls.sim[i,] <- nnls.sim[i,]/sum(nnls.sim[i,])
}
colnames(nnls.sim) <- colnames(X)

betas.sim <- as.matrix(read.csv('simulation_betas.csv',row.names=1))
groups.sim <- read.csv('simulation_groups.csv',row.names=1,stringsAsFactors=F)
groups.sim[groups.sim=='Left atrium'] <- 'Heart'

nnls.sim <- as.matrix(read.csv('simulation_nnls.csv',row.names=1,check.names=F))

rownames(nnls.sim) <- rownames(betas.sim)
betas.sim.all <- matrix(0,nrow(nnls.sim),ncol(nnls.sim))
colnames(betas.sim.all) <- colnames(nnls.sim)
rownames(betas.sim.all) <- rownames(nnls.sim)

write.csv(nnls.sim,'simulation_nnls.csv',row.names=T)

for (i in 1:nrow(betas.sim.all)){
	betas.sim.all[i,match(groups.sim[i,],colnames(betas.sim.all))] <- betas.sim[i,]
}

cors <- sapply(1:nrow(nnls.sim),function(i) cor(betas.sim.all[i,],nnls.sim[i,]))

betas.diff <- nnls.sim-betas.sim.all
betas.diff.abs <- abs(betas.diff) 
betas.diff.abs.mean <- apply(betas.diff.abs,2,mean)

ord <- order(betas.diff.abs.mean,decreasing=T)
betas.diff.v <- as.vector(betas.diff[,ord])
betas.diff.abs.v <- as.vector(betas.diff.abs[,ord])

x <- rep((1:ncol(betas.diff)),each=nrow(betas.diff))
lab <- rep(1:nrow(betas.diff),times=ncol(betas.diff))

pdf('simulation_test_confusion_plot.pdf')
par(mar=c(15, 4, 4, 2) + 0.1)
stripchart(as.data.frame(betas.diff[,ord]),las=2,pch=NA,vertical=T)
text(x[betas.diff.abs.v>0.05],betas.diff.v[betas.diff.abs.v>0.05],lab[betas.diff.abs.v>0.05])
dev.off()


#stripchart(as.data.frame(betas.diff.abs[,order(betas.diff.abs.mean,decreasing=T)]),las=2,pch=1,vertical=T)
stripchart(as.data.frame(betas.diff[,order(betas.diff.abs.mean,decreasing=T)]),las=2,pch=1,vertical=T)

pdf('simulation_test.pdf')
for (i in 1:nrow(betas.sim.all)){
	plot(betas.sim.all[i,],nnls.sim[i,],type='n',xlab='Simulated betas',ylab='Predicted betas',main=paste0('Set',i))
	text(betas.sim.all[i,],nnls.sim[i,],labels=colnames(betas.sim.all))
}
dev.off()

# look at mistakes from zero
false.pos <- betas.sim.all==0 & nnls.sim>0.01
false.pos.group <- apply(false.pos,2,sum)

# find anti-correlated tissues
library(Hmisc)
cor <- rcorr(betas.diff,type='pearson')
library(psych)
cor.psych <- corr.test(betas.diff, use = "pairwise",method="pearson",adjust="fdr",alpha=.05)
cor.T <- (cor.psych$p < 0.05 & cor.psych$r<0)*1
for (i in 1:nrow(cor.T)){
	for (j in 1:i){
		cor.T[i,j] <- cor.T[j,i]
	}
}
library(gplots)
pdf('anti-correlated_residuals.pdf')
heatmap.2(cor.T,trace='none',main='anti-correlated residuals',margins=c(10,10))
dev.off()
