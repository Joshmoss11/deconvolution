source('nnls_plots_functions.R')
deconv.dir <- '/mnt/lustre/hms-01/fs01/joshua.moss/dor/deconvolution'
nnls.dir <- file.path(deconv.dir,'nnls_select')
plots.dir <- file.path(deconv.dir,'plots_select')
dir.create(plots.dir,showWarnings=F)
matrix.f <- file.path(nnls.dir,'nnls_betas_norm.csv')
pies.f <- file.path(plots.dir,'betas_pies.pdf')
#heatmap.f <- file.path(plots.dir,'betas_heatmap.pdf')

matrix.f.w.c <- file.path(nnls.dir,'nnls_betas_norm_wi_cancer.csv')
pies.f.w.c <- file.path(plots.dir,'betas_pies_wi_cancer.pdf')
#heatmap.f.w.c <- file.path(plots.dir,'betas_heatmap_wi_cancer.pdf')

make_nnls_pies(matrix.f,pies.f)
#make_nnls_heatmap(matrix.f,heatmap.f)

make_nnls_pies(matrix.f.w.c,pies.f.w.c)
#make_nnls_heatmap(matrix.f.w.c,heatmap.f.w.c)

matrix.f.ve <- file.path(nnls.dir,'nnls_betas_norm_wi_vasc_end.csv')
pies.f.ve <- file.path(plots.dir,'betas_pies_wi_vasc_end.pdf')

matrix.f.w.c.ve <- file.path(nnls.dir,'nnls_betas_norm_wi_cancer_wi_vasc_end.csv')
pies.f.w.c.ve <- file.path(plots.dir,'betas_pies_wi_cancer_wi_vasc_end.pdf')

make_nnls_pies(matrix.f.ve,pies.f.ve)
make_nnls_pies(matrix.f.w.c.ve,pies.f.w.c.ve)
