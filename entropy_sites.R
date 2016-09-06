entropy_sites <- function(X,num_sites){

    r <- dim(X)[1]; c<- dim(X)[2]
    X.norm <- matrix(NA,r,c)
    sites.idx.1 <- matrix(NA,num_sites,c)
    sites.idx.2 <- matrix(NA,num_sites,c)
    # normalize each cpg
    for (i in 1:r){
        X.norm[i,] <- X[i,]/sum(X[i,])
    }
    # get top for each tissue
    for (j in 1:c){
        sites.idx.1[,j] <- order(X.norm[,j],decreasing = T,na.last = T)[1:num_sites]
    }
    ## and in other direction
    X.rev <- 1-X
    # normalize each cpg
    for (i in 1:r){
        X.norm[i,] <- X.rev[i,]/sum(X.rev[i,])
    }
    # get top for each tissue
    for (j in 1:c){
      sites.idx.2[,j] <- order(X.norm[,j],decreasing = T,na.last = T)[1:num_sites]
    }

    sites.vector <- c(as.vector(sites.idx.1),as.vector(sites.idx.2))
    return(sites.vector)
}
