nrows=-1
matrices.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution/matrices'
data.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data'
betas.dir <- file.path(data.dir,'lab_data/betas/by_tissue')

files.to.read <- c('cfDNA.csv','Mix.csv','Leukocyte.csv','Leukocytes_cord_blood.csv')
for (f in files.to.read){
	print(paste0('Reading ',f))
	d <- read.csv(file.path(betas.dir,f),header=T,row.names=1,check.names=F,nrow=nrows)
	if (f == files.to.read[1]){
		d.all <- d
	} else {
		d.all <- merge(d.all,d,by=0);rownames(d.all) <- d.all[,1]; d.all <- d.all[,2:ncol(d.all)]
	}
}

# change names to ids
annot.f <- file.path(data.dir,'lab_data','sample_sheet.csv')
annot <- read.csv(annot.f,stringsAsFactors=F,check.names=F)

new.names <- annot$Sample[match(colnames(d.all),annot$barcode)]
colnames(d.all) <- new.names

#write to file
print(paste0('Writing test matrix to file'))
write.csv(d.all,file.path(matrices.dir,'test_samples.csv'),row.names=T,quote=F)
