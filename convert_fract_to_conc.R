convert_fract_to_conc <- function(fract.f,conc.f,cfdna.conc.f = '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution/annotations/cfdna_conc.csv'){
	fract <- as.matrix(read.csv(fract.f,row.names=1,check.names=F))
	cfdna.conc <- read.csv(cfdna.conc.f)
	cfdna.conc <- cfdna.conc[match(rownames(fract),cfdna.conc[,1]),2]
	conc <- apply(fract,2,function(x) x*cfdna.conc)
	write.csv(conc,conc.f,row.names=T,quote=F)
}

deconv.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution/'
nnls.dir <- file.path(deconv.dir,'nnls')
conc.dir <- file.path(deconv.dir,'conc')
dir.create(conc.dir,showWarnings=F)
fract.f <- file.path(nnls.dir,'nnls_betas_norm_wi_vasc_end.csv')
conc.f <- file.path(conc.dir,'nnls_betas_norm_wi_vasc_end.csv')
convert_fract_to_conc(fract.f,conc.f)
