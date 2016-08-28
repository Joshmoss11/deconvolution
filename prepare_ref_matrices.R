nrows=-1
matrices.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution/matrices'
data.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data'

#create ref matrix

#without cancer

print(paste0('Reading in lab data'))
lab_data.bad_cols <- c('Leukocyte','Leukocytes_cord_blood','Mix','Negative','cfDNA','Fat','Liver')
lab_data.means.f <- file.path(data.dir,'lab_data/summary/by_statistics/means.csv')
lab_data.means <- read.csv(lab_data.means.f,header=T,row.names=1,check.names=F,nrow=nrows)
lab_data.means <- lab_data.means[,!(colnames(lab_data.means) %in% lab_data.bad_cols)]

print(paste0('Reading in GEO'))
geo.bad_cols <- c('Cord_blood','Endothelium','Granulocytes','Pancreas','PBMC','Sc_fat','Whole_blood','Whole_maternal_blood') 
geo.means.f <- file.path(data.dir,'GEO/summary/by_statistics/means.csv')
geo.means <- read.csv(geo.means.f,header=T,row.names=1,check.names=F,nrow=nrows)
geo.means <- geo.means[,!(colnames(geo.means) %in% geo.bad_cols)]

print(paste0('Reading in TCGA'))
tcga.bad_cols <- c('Kidney_2','Liver','Lung_2','Pancreas','Brain')
tcga.means.f <- file.path(data.dir,'TCGA/summary/by_statistics/Solid_Tissue_Normal/means_Solid_Tissue_Normal.csv')
tcga.means <- read.csv(tcga.means.f,header=T,row.names=1,check.names=F,nrow=nrows)
colnames(tcga.means) <- c('Bladder','Breast','Cervix','Colon','Esophagus','Brain','Head_and_kneck','Kidney','Kidney_2','Liver',
							'Lung','Lung_2','Pancreas','Testes','Prostate','Rectum','Mesenchyme','Skin','Stomach','Thyroid','Thymus','Uterus')
tcga.means <- tcga.means[,!(colnames(tcga.means) %in% tcga.bad_cols)]

ref <- merge(lab_data.means,geo.means,by=0); rownames(ref) <- ref[,1]; ref <- ref[,2:ncol(ref)]
ref <- merge(ref,tcga.means,by=0); rownames(ref) <- ref[,1]; ref <- ref[,2:ncol(ref)]

print(paste0('Writing reference matrix to file'))
write.csv(ref,file.path(matrices.dir,'ref_no_cancer.csv'),row.names=T,quote=F)

#with cancer
print(paste0('Reading in cancer data'))
tcga.cancer.means.f <- file.path(data.dir,'TCGA/summary/by_statistics/Primary_Tumor/means_Primary_Tumor.csv')
tcga.cancer.means <- read.csv(tcga.cancer.means.f,header=T,row.names=1,check.names=F,nrow=nrows)
print(paste0('Writing reference matrix (wih cancer) to file'))
ref.w.cancer <- merge(ref,tcga.cancer.means,by=0); rownames(ref.w.cancer) <- ref.w.cancer[,1]; ref.w.cancer <- ref.w.cancer[,2:ncol(ref.w.cancer)]
write.csv(ref.w.cancer,file.path(matrices.dir,'ref_wi_cancer.csv'),row.names=T,quote=F)
