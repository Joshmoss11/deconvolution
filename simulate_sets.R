convert_tcga_names <- function(tcga_names_query){
    tcga.names <- c('BLCA',  'LUSC','BRCA',  'PAAD','CESC',  'PCPG','COAD',  'PRAD',
		'ESCA',  'READ','GBM',   'SARC','HNSC',  'SKCM','KIRC',  'STAD','KIRP',  'THCA','LIHC',  'THYM','LUAD',  'UCEC')
    tcga.names.tissue <- c('Bladder','Lung_2','Breast','Pancreas','Cervix','Testes','Colon','Prostate',
		'Esophagus','Rectum','Brain','Mesenchyme','Head_and_kneck','Skin','Kidney','Stomach','Kidney_2','Thyroid','Liver','Thymus','Lung','Uterus')
    return(tcga.names.tissue[match(tcga_names_query,tcga.names)])
}
get_num_samples <- function(csv){
    x <- read.csv(csv,row.names=1,nrows=10)
    return(ncol(x))
}
read_random_sample <- function(csv){
    s <- get_num_samples(csv)
    cl <- c('character',rep('NULL',s))
    cl[sample(1:s,1)+1] <- 'numeric'
    x <- read.csv(csv,row.names=1,colClasses=cl)
    return(x)
}
simulate_set <- function(me.groups,me.files,num.groups){
    groups.idx <- sample(1:length(me.groups),num.groups)
    set.groups <- me.groups[groups.idx]
    set.files <- me.files[groups.idx]
    set.betas <- runif(num.groups,0,1)
    set.betas <- set.betas/sum(set.betas)
    for (j in 1:num.groups){
        cur.y <- read_random_sample(set.files[j])
        if (j==1) {y <- cur.y} else {
            y <- merge(y,cur.y,by=0); rownames(y) <- y[,1]; y <- y[,2:ncol(y)]
        }
    }
    y <- as.matrix(y) %*% set.betas
    return(list(y,set.groups,set.betas))
}
data.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data'
lab_data.dir <- file.path(data.dir,'lab_data/betas/by_tissue')
geo.dir <- file.path(data.dir,'GEO/betas/by_tissue')
tcga.dir <- file.path(data.dir,'TCGA/betas/by_status')
tcga.healthy.dir <- file.path(tcga.dir,'Solid_Tissue_Normal')
tcga.cancer.dir <- file.path(tcga.dir,'Primary_Tumor')
sim.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution/simulations'
dir.create(sim.dir,showWarnings=F)
lab_data.files <- list.files(lab_data.dir)[grep('.csv',list.files(lab_data.dir))]
lab_data.tissues <- substr(lab_data.files,1,nchar(lab_data.files)-nchar('.csv'))
lab_data.tissues.bad <- c('Leukocyte','Leukocytes_cord_blood','Mix','Negative','cfDNA','Fat','Liver')
lab_data.files <- file.path(lab_data.dir,lab_data.files[!(lab_data.tissues %in% lab_data.tissues.bad)])
lab_data.tissues <- lab_data.tissues[!(lab_data.tissues %in% lab_data.tissues.bad)]

geo.files <- list.files(geo.dir)[grep('.csv',list.files(geo.dir))]
geo.tissues <- substr(geo.files,1,nchar(geo.files)-nchar('.csv'))
geo.tissues.bad <- c('Granulocytes','Liver','Omentum','PBMC','Pancreas','Sc_fat','Whole_blood','Whole_maternal_blood')
geo.files <- file.path(geo.dir,geo.files[!(geo.tissues %in% geo.tissues.bad)])
geo.tissues <- geo.tissues[!(geo.tissues %in% geo.tissues.bad)]

tcga.healthy.tissues.bad <- c('Thymus','Kidney_2','Liver','Lung_2','Pancreas','Brain','Skin')

tcga.healthy.files <- list.files(tcga.healthy.dir)[grep('.csv',list.files(tcga.healthy.dir))]
tcga.healthy.tissues <- convert_tcga_names(substr(tcga.healthy.files,1,nchar(tcga.healthy.files)-nchar('_Solid_Tissue_Normal.csv')))
tcga.healthy.files <- file.path(tcga.healthy.dir,tcga.healthy.files[!(tcga.healthy.tissues %in% tcga.healthy.tissues.bad)])
tcga.healthy.tissues <- tcga.healthy.tissues[!(tcga.healthy.tissues %in% tcga.healthy.tissues.bad)]

tcga.cancer.files <- list.files(tcga.cancer.dir)[grep('.csv',list.files(tcga.cancer.dir))]
tcga.cancer.tissues <- substr(tcga.cancer.files,1,nchar(tcga.cancer.files)-nchar('_Primary_Tumor.csv'))
tcga.cancer.files <- file.path(tcga.cancer.dir,tcga.cancer.files)

healthy.files <- c(lab_data.files,geo.files,tcga.healthy.files)
all.files <- c(healthy.files,tcga.cancer.files)

healthy.tissues <- c(lab_data.tissues,geo.tissues,tcga.healthy.tissues)
all.tissues <- c(healthy.tissues,tcga.cancer.tissues)

num.groups <- 10
num.sets <- 1000

sets.groups <- matrix(NA,num.sets,num.groups)
sets.betas <- matrix(NA,num.sets,num.groups)

for (i in 1:num.sets) {
	print(paste0('Simulating dataset # ',i,'/',num.sets))
	set <- simulate_set(all.tissues,all.files,num.groups)
	y <- set[[1]]
	sets.groups[i,] <- set[[2]]
	sets.betas[i,] <- set[[3]]	
	if (i==1) {sets.y <- data.frame(y)} else {
		sets.y <- merge(sets.y,y,by=0); rownames(sets.y) <- sets.y[,1]; sets.y <- sets.y[,2:ncol(sets.y)]
	}
}
colnames(sets.y) <- paste0('set',seq(1,num.sets,1))
colnames(sets.groups) <- paste0('group',seq(1,num.groups,1))
colnames(sets.betas) <- paste0('group',seq(1,num.groups,1))
rownames(sets.groups) <- paste0('set',seq(1,num.sets,1))
rownames(sets.betas) <- paste0('set',seq(1,num.sets,1))

write.csv(sets.groups,file.path(sim.dir,'simulation_groups.csv'),quote=F,row.names=T)
write.csv(sets.betas,file.path(sim.dir,'simulation_betas.csv'),quote=F,row.names=T)
write.csv(sets.y,file.path(sim.dir,'simulation_y.csv'),quote=F,row.names=T)
