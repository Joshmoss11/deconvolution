source('nnls_matrix_functions.R')
deconv.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution'
matrices.dir <- file.path(deconv.dir,'matrices')
nnls.dir <- file.path(deconv.dir,'nnls_select')
dir.create(nnls.dir,showWarnings=F)
X.f <- file.path(matrices.dir,'ref_no_cancer.csv')
Y.f <- file.path(matrices.dir,'test_samples.csv')
out.f <- file.path(nnls.dir,'nnls_betas.csv')
out.norm.f <- file.path(nnls.dir,'nnls_betas_norm.csv')
X.f.w.c <- file.path(matrices.dir,'ref_wi_cancer.csv')
out.f.w.c <- file.path(nnls.dir,'nnls_betas_wi_cancer.csv')
out.norm.f.w.c <- file.path(nnls.dir,'nnls_betas_norm_wi_cancer.csv')

out.f.ve <- file.path(nnls.dir,'nnls_betas_wi_vasc_end.csv')
out.norm.f.ve <- file.path(nnls.dir,'nnls_betas_norm_wi_vasc_end.csv')
out.f.w.c.ve <- file.path(nnls.dir,'nnls_betas_wi_cancer_wi_vasc_end.csv')
out.norm.f.w.c.ve <- file.path(nnls.dir,'nnls_betas_norm_wi_cancer_wi_vasc_end.csv')

run_nnls_matrix(X.f,Y.f,out.f,out.norm.f,entropy_sites=T)
run_nnls_matrix(X.f.w.c,Y.f,out.f.w.c,out.norm.f.w.c,entropy_sites=T)

run_nnls_matrix(X.f,Y.f,out.f.ve,out.norm.f.ve,rm.cols=c('Skin','Liver','Omentum'),entropy_sites=T)
run_nnls_matrix(X.f.w.c,Y.f,out.f.w.c.ve,out.norm.f.w.c.ve,rm.cols=c('Skin','Liver','Omentum'),entropy_sites=T)
