#######################################################################################
#                            Simulation to test Sparse PCA                            #
# Test accuracy of sparse covariance matrix estimation
# Case 1: dense
# Case 2: sparse
#######################################################################################

################                 Requare packages                      ################
library(ggplot2)
library(PMA)
#library(rospca)
library(pipeR)
library(pracma) # to generate orthonormal matices
library(Matrix)
################                Sample generation function             ################

generateX <- function(eigenvec, eigenval, n) {
        p = length(eigenval)
        Xmat <- matrix(0, n, p)
        for (i in 1:n){
                pcascore <- rnorm(p, mean = 0, sd = sqrt(eigenval))
                Xmat[i, ] <- eigenvec%*%pcascore # Xi
        }
        return(Xmat) # return n by p matrix
}

################         Function to get prop.var.explained            ################

prop.var.ex <- function(Xmatrix, K, method) {
        # Xmatrix = sample matirx
        # K = number of factors in PMD return
        # method = "sparse" or "ordinary"
        if (method == "sparse") {
                sumab.cv <- SPC.cv(Xmatrix, sumabsvs = seq(1.2, sqrt(50), length.out = 15))$bestsumabsv
                prop.var <- SPC(Xmatrix, sumabsv = sumab.cv, K= K)$prop.var.explained
        } else {
                eigenvalue <- eigen(t(Xmatrix)%*%Xmatrix/ncol(Xmatrix))$values
                prop.var <- cumsum(eigenvalue[1:K])/sum(eigenvalue)
        }
        return(prop.var)
}

################    Function to get prop.var.explained np times        ################

prop.matrix <- function (np, K, method, eigenvec, eigenval, n) {
        # np = number of repetitions
        res.matrix <- matrix(0, np, K)
        for (j in 1:np) {
                Xmat <- generateX(eigenvec = eigenvec, eigenval = eigenval, n = n)
                res.matrix[j, ] <- prop.var.ex(Xmatrix = Xmat, K = K, method = method)
        }
        return(res.matrix)
}
#######################################################################################
################                CASE 1: dense eigenvectors              ###############
#######################################################################################
# number of variables, p = 50
# sample size, n = 20, 40, 80, 100, 200
# number of repetitions for each case: 20

################                 Simulation Setting                    ################
set.seed(20180422)
EigenVect1 <- randortho(50, type = "orthonormal")
EigenValue1 <- c(seq(5,0.5,length.out = 10),rep(0.1, 40))
True.Var.Prop <- cumsum(EigenValue1[1:15])/sum(EigenValue1)

################                    Sparse PCA method                  ################

prop.var.dense.20 <- prop.matrix(np = 20, K = 15, method = "sparse", eigenvec = EigenVect1, 
                                 eigenval = EigenValue1, n = 20) %>>% colMeans()
prop.var.dense.40 <- prop.matrix(np = 20, K = 15, method = "sparse", eigenvec = EigenVect1,
                                  eigenval = EigenValue1, n = 40) %>>% colMeans()
prop.var.dense.80 <- prop.matrix(np = 20, K = 15, method = "sparse", eigenvec = EigenVect1,
                                 eigenval = EigenValue1, n = 80) %>>% colMeans()
prop.var.dense.100 <- prop.matrix(np = 20, K = 15, method = "sparse", eigenvec = EigenVect1,
                                  eigenval = EigenValue1, n = 100) %>>% colMeans()
prop.var.dense.10 <- prop.matrix(np = 20, K = 10, method = "sparse", eigenvec = EigenVect1,
                                 eigenval = EigenValue1, n = 10) %>>% colMeans()

###############            Eigen decomposition method                  ################

prop.var.dense.pca.20 <- prop.matrix(np = 20, K = 15, method = "ordinary", eigenvec = EigenVect1,
                                  eigenval = EigenValue1, n = 20) %>>% colMeans()
prop.var.dense.pca.40 <- prop.matrix(np = 20, K = 15, method = "ordinary", eigenvec = EigenVect1,
                                  eigenval = EigenValue1, n = 40) %>>% colMeans()
prop.var.dense.pca.80 <- prop.matrix(np = 20, K = 15, method = "ordinary", eigenvec = EigenVect1,
                                  eigenval = EigenValue1, n = 80) %>>% colMeans()
prop.var.dense.pca.100 <- prop.matrix(np = 20, K = 15, method = "ordinary", eigenvec = EigenVect1,
                                  eigenval = EigenValue1, n = 100) %>>% colMeans()
prop.var.dense.pca.10 <- prop.matrix(np = 20, K = 15, method = "ordinary", eigenvec = EigenVect1,
                                  eigenval = EigenValue1, n = 10) %>>% colMeans()

################             Plot the prop.var.explained               ################
DenseEigen.SparsePCA <- data.frame(prop.var.exp = c(True.Var.Prop, prop.var.dense.20,prop.var.dense.40,prop.var.dense.80,prop.var.dense.100,prop.var.dense.10,NA,NA,NA,NA,NA),
                                   K = rep (1:15, 6),
                                   SampleSize = rep(c("True.Var.Prop", "n = 20" , "n = 40" , "n = 80" , "n = 100" , "n = 10"), each = 15))
par(mar=c(2,2,2,2))
ggplot(DenseEigen.SparsePCA,aes(x = K,y = prop.var.exp,colour = SampleSize)) +
        geom_point() +
        geom_line() +
        ggtitle("Dense Egienvector solved by Spare PCA method")+
        labs(x = "Number of Factors", y = "Proportion of Variance Explained")

DenseEigen.StandPCA <- data.frame(prop.var.exp = c(True.Var.Prop, prop.var.dense.pca.20,prop.var.dense.pca.40,prop.var.dense.pca.80,prop.var.dense.pca.100,prop.var.dense.pca.10),
                                   K = rep (1:15, 6),
                                   SampleSize = rep(c("True.Var.Prop", "n = 20" , "n = 40" , "n = 80" , "n = 100" , "n = 10"), each = 15))
ggplot(DenseEigen.StandPCA,aes(x = K,y = prop.var.exp,colour = SampleSize)) +
        geom_point() +
        geom_line() +
        ggtitle("Dense Egienvector solved by regular PCA method")+
        labs(x = "Number of Factors", y = "Proportion of Variance Explained")


#######################################################################################
################                CASE 2: sparse eigenvectors              ##############
#######################################################################################
# number of variables, p = 50
# sample size, n = 20, 40, 80, 100, 200
# number of repetitions for each case: 20
set.seed(20180421)
EigenVect <- randortho(50, type = "orthonormal")
EigenValue2 <- c(seq(5,0.5,length.out = 10),rep(0.1, 40))

Xmat <- generateX(eigenvec = EigenVect, eigenval = EigenValue, n = 200)
res <- SPC(Xmat, sumabsv = 2, K= 50, orth =T)
EigenVect2 <- res$v
NumNonZero <- c(6, 6, 7, 9, 7, 6, 8, 11, 9, 7, 8, 7, 9, 15, 7, 12, 6, 7, 6, 8, 8, 
                9, 8, 8, 10, 6, 7, 8, 9, 9, 7, 10, 10, 9, 11, 11, 12, 14, 6, 9, 
                14, 8, 13, 8, 6, 6, 12, 6, 6, 11)

################                    Sparse PCA method                  ################

prop.var.sparse.20 <- prop.matrix(np = 20, K = 15, method = "sparse", eigenvec = EigenVect2, 
                                 eigenval = EigenValue2, n = 20) %>>% colMeans()
prop.var.sparse.40 <- prop.matrix(np = 20, K = 15, method = "sparse", eigenvec = EigenVect2,
                                 eigenval = EigenValue2, n = 40) %>>% colMeans()
prop.var.sparse.80 <- prop.matrix(np = 20, K = 15, method = "sparse", eigenvec = EigenVect2,
                                 eigenval = EigenValue2, n = 80) %>>% colMeans()
prop.var.sparse.100 <- prop.matrix(np = 20, K = 15, method = "sparse", eigenvec = EigenVect2,
                                  eigenval = EigenValue2, n = 100) %>>% colMeans()
prop.var.sparse.10 <- prop.matrix(np = 20, K = 10, method = "sparse", eigenvec = EigenVect2,
                                 eigenval = EigenValue2, n = 10) %>>% colMeans()

###############            Eigen decomposition method                  ################

prop.var.sparse.pca.20 <- prop.matrix(np = 20, K = 15, method = "ordinary", eigenvec = EigenVect2,
                                     eigenval = EigenValue2, n = 20) %>>% colMeans()
prop.var.sparse.pca.40 <- prop.matrix(np = 20, K = 15, method = "ordinary", eigenvec = EigenVect2,
                                     eigenval = EigenValue2, n = 40) %>>% colMeans()
prop.var.sparse.pca.80 <- prop.matrix(np = 20, K = 15, method = "ordinary", eigenvec = EigenVect2,
                                     eigenval = EigenValue2, n = 80) %>>% colMeans()
prop.var.sparse.pca.100 <- prop.matrix(np = 20, K = 15, method = "ordinary", eigenvec = EigenVect2,
                                      eigenval = EigenValue2, n = 100) %>>% colMeans()
prop.var.sparse.pca.10 <- prop.matrix(np = 20, K = 15, method = "ordinary", eigenvec = EigenVect2,
                                     eigenval = EigenValue2, n = 10) %>>% colMeans()

################             Plot the prop.var.explained               ################
SparseEigen.SparsePCA <- data.frame(prop.var.exp = c(True.Var.Prop, prop.var.sparse.20,prop.var.sparse.40,prop.var.sparse.80,prop.var.sparse.100,prop.var.sparse.10,NA,NA,NA,NA,NA),
                                   K = rep (1:15, 6),
                                   SampleSize = rep(c("True.Var.Prop", "n = 20" , "n = 40" , "n = 80" , "n = 100" , "n = 10"), each = 15))
par(mar=c(2,2,2,2))
ggplot(SparseEigen.SparsePCA,aes(x = K,y = prop.var.exp,colour = SampleSize)) +
        geom_point() +
        geom_line() +
        ggtitle("Sparse Egienvector solved by Spare PCA method")+
        labs(x = "Number of Factors", y = "Proportion of Variance Explained")

SparseEigen.StandPCA <- data.frame(prop.var.exp = c(True.Var.Prop, prop.var.sparse.pca.20,prop.var.sparse.pca.40,prop.var.sparse.pca.80,prop.var.sparse.pca.100,prop.var.sparse.pca.10),
                                  K = rep (1:15, 6),
                                  SampleSize = rep(c("True.Var.Prop", "n = 20" , "n = 40" , "n = 80" , "n = 100" , "n = 10"), each = 15))
ggplot(SparseEigen.StandPCA,aes(x = K,y = prop.var.exp,colour = SampleSize)) +
        geom_point() +
        geom_line() +
        ggtitle("Sparse Egienvector solved by regular PCA method")+
        labs(x = "Number of Factors", y = "Proportion of Variance Explained")



#######################################################################################
##############       Estimation Error -- CASE 1: dense eigenvectors      ##############
#######################################################################################

################         Function to get estimation error              ################

prop.error <- function(Xmatrix, K, method, trueSigma) {
        # Xmatrix = sample matirx
        # K = number of factors in PMD return
        # method = "sparse" or "ordinary"
        if (method == "sparse") {
                sumab.cv <- SPC.cv(Xmatrix, sumabsvs = seq(1.2, sqrt(50), length.out = 15))$bestsumabsv
                res <- SPC(Xmatrix, sumabsv = sumab.cv, K= K)
                sigma.appr <- (res$v)%*%diag(res$d)%*%t(res$u)%*%(res$u)%*%(diag(res$d))%*%t(res$v)/nrow(Xmatrix)
        } else {
                out <- eigen((t(Xmatrix)%*%Xmatrix)/nrow(Xmatrix))
                sigma.appr <- (out$vectors[, 1:K])%*%diag(out$values[1:K])%*%t(out$vectors[, 1:K])
        }
        prop <- norm(trueSigma - sigma.appr, "F")/norm(trueSigma, "F")
}

################    Function to get prop.var.explained np times        ################

error.vec <- function (np, K, method, eigenvec, eigenval, n, trueSigma) {
        # np = number of repetitions
        res.vec <- rep(0, np)
        for (j in 1:np) {
                Xmat <- generateX(eigenvec = eigenvec, eigenval = eigenval, n = n)
                res.vec[j]<- prop.error(Xmatrix = Xmat, K = K, method = method, trueSigma = trueSigma)
        }
        return(res.vec)
}

trueSigma1 <- EigenVect1%*%diag(EigenValue1)%*%t(EigenVect1)

True.Error.Prop <- norm(trueSigma1 - EigenVect1[,1:10]%*%diag(EigenValue1[1:10])%*%t(EigenVect1[,1:10]), "F")/norm(trueSigma1, "F")

err.dense.10 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect1,
                          eigenval = EigenValue1, n = 10, trueSigma = trueSigma1)
err.dense.20 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect1,
                          eigenval = EigenValue1, n = 20, trueSigma = trueSigma1)
err.dense.40 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect1,
                          eigenval = EigenValue1, n = 40, trueSigma = trueSigma1)
err.dense.80 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect1,
                          eigenval = EigenValue1, n = 80, trueSigma = trueSigma1)
err.dense.100 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect1,
                          eigenval = EigenValue1, n = 100, trueSigma = trueSigma1)
err.dense.200 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect1,
                           eigenval = EigenValue1, n = 200, trueSigma = trueSigma1)
err.dense.400 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect1,
                           eigenval = EigenValue1, n = 400, trueSigma = trueSigma1)

err.dense.pca.10 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect1,
                          eigenval = EigenValue1, n = 10, trueSigma = trueSigma1)
err.dense.pca.20 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect1,
                          eigenval = EigenValue1, n = 20, trueSigma = trueSigma1)
err.dense.pca.40 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect1,
                          eigenval = EigenValue1, n = 40, trueSigma = trueSigma1)
err.dense.pca.80 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect1,
                          eigenval = EigenValue1, n = 80, trueSigma = trueSigma1)
err.dense.pca.100 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect1,
                           eigenval = EigenValue1, n = 100, trueSigma = trueSigma1)
err.dense.pca.200 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect1,
                           eigenval = EigenValue1, n = 200, trueSigma = trueSigma1)
err.dense.pca.400 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect1,
                               eigenval = EigenValue1, n = 400, trueSigma = trueSigma1)



#######################################################################################
##############       Estimation Error -- CASE 2: sparse eigenvectors     ##############
#######################################################################################
trueSigma2 <- EigenVect2%*%diag(EigenValue2)%*%t(EigenVect2)

True.Error.Prop <- norm(trueSigma2 - EigenVect2[,1:10]%*%diag(EigenValue2[1:10])%*%t(EigenVect2[,1:10]), "F")/norm(trueSigma2, "F")

err.sparse.10 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect2,
                          eigenval = EigenValue2, n = 10, trueSigma = trueSigma2)
err.sparse.20 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect2,
                          eigenval = EigenValue2, n = 20, trueSigma = trueSigma2)
err.sparse.40 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect2,
                          eigenval = EigenValue2, n = 40, trueSigma = trueSigma2)
err.sparse.80 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect2,
                          eigenval = EigenValue2, n = 80, trueSigma = trueSigma2)
err.sparse.100 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect2,
                           eigenval = EigenValue2, n = 100, trueSigma = trueSigma2)
err.sparse.200 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect2,
                           eigenval = EigenValue2, n = 200, trueSigma = trueSigma2)
err.sparse.400 <- error.vec(np = 20, K = 10, method = "sparse", eigenvec = EigenVect2,
                            eigenval = EigenValue2, n = 400, trueSigma = trueSigma2)




err.sparse.pca.10 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect2,
                          eigenval = EigenValue2, n = 10, trueSigma = trueSigma2)
err.sparse.pca.20 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect2,
                          eigenval = EigenValue2, n = 20, trueSigma = trueSigma2)
err.sparse.pca.40 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect2,
                          eigenval = EigenValue2, n = 40, trueSigma = trueSigma2)
err.sparse.pca.80 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect2,
                          eigenval = EigenValue2, n = 80, trueSigma = trueSigma2)
err.sparse.pca.100 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect2,
                           eigenval = EigenValue2, n = 100, trueSigma = trueSigma2)
err.sparse.pca.200 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect2,
                           eigenval = EigenValue2, n = 200, trueSigma = trueSigma2)
err.sparse.pca.400 <- error.vec(np = 20, K = 10, method = "ordinary", eigenvec = EigenVect2,
                                eigenval = EigenValue2, n = 400, trueSigma = trueSigma2)

SparseEigen.SparsePCA.err <- rbind(err.sparse.10,err.sparse.20,err.sparse.40,err.sparse.80,
                                   err.sparse.100,err.sparse.200,err.sparse.400)
SparseEigen.SparsePCA.err.mean <- apply(SparseEigen.SparsePCA.err, 1, mean)
SparseEigen.SparsePCA.err.sd <- apply(SparseEigen.SparsePCA.err, 1, sd)


SparseEigen.reguPCA.err <- rbind(err.sparse.pca.10,err.sparse.pca.20,err.sparse.pca.40,err.sparse.pca.80,
                                   err.sparse.pca.100,err.sparse.pca.200,err.sparse.pca.400)
SparseEigen.reguPCA.err.mean <- apply(SparseEigen.reguPCA.err, 1, mean)
SparseEigen.reguPCA.err.sd <- apply(SparseEigen.reguPCA.err, 1, sd)

denseEigen.SparsePCA.err <- rbind(err.dense.10,err.dense.20,err.dense.40,err.dense.80,
                                   err.dense.100,err.dense.200,err.dense.400)
denseEigen.SparsePCA.err.mean <- apply(denseEigen.SparsePCA.err, 1, mean)
denseEigen.SparsePCA.err.sd <- apply(denseEigen.SparsePCA.err, 1, sd)


denseEigen.reguPCA.err <- rbind(err.dense.pca.10,err.dense.pca.20,err.dense.pca.40,err.dense.pca.80,
                                  err.dense.pca.100,err.dense.pca.200,err.dense.pca.400)
denseEigen.reguPCA.err.mean <- apply(denseEigen.reguPCA.err, 1, mean)
denseEigen.reguPCA.err.sd <- apply(denseEigen.reguPCA.err, 1, sd)


SparseEigen.SparsePCA.err.df<-data.frame(SampleSize = rep(c(10, 20, 40, 80, 100, 200, 400),2),
                                         Err.percen = c(SparseEigen.SparsePCA.err.mean, SparseEigen.reguPCA.err.mean),
                                         Method = rep(c("Sparse PCA", "Regular PCA"), each = 7),
                                         CI.upper1 = rep(SparseEigen.SparsePCA.err.mean + 2*SparseEigen.SparsePCA.err.sd, 2),
                                         CI.lower1 = rep(SparseEigen.SparsePCA.err.mean - 2*SparseEigen.SparsePCA.err.sd, 2),
                                         CI.upper2 = rep(SparseEigen.reguPCA.err.mean + 2*SparseEigen.reguPCA.err.sd, 2),
                                         CI.lower2 = rep(SparseEigen.reguPCA.err.mean - 2*SparseEigen.reguPCA.err.sd, 2))
                                         
                                        
qplot(SampleSize, Err.percen, data=SparseEigen.SparsePCA.err.df, geom="line", col = Method,
      main="Error Percentage in Case 2-Sparse Eigenvector",size=I(2))+
        geom_ribbon(aes(ymin=CI.lower1,ymax=CI.upper1), alpha=0.1,col = "darkgreen")+
        geom_ribbon(aes(ymin=CI.lower2,ymax=CI.upper2), alpha=0.1,col = "red")+
        geom_hline(yintercept=True.Error.Prop,linetype="dashed", color = "black",size=1)


        
denseEigen.SparsePCA.err.df<-data.frame(SampleSize = rep(c(10, 20, 40, 80, 100, 200, 400),2),
                                         Err.percen = c(denseEigen.SparsePCA.err.mean, denseEigen.reguPCA.err.mean),
                                         Method = rep(c("Sparse PCA", "Regular PCA"), each = 7),
                                         CI.upper1 = rep(denseEigen.SparsePCA.err.mean + 2*denseEigen.SparsePCA.err.sd, 2),
                                         CI.lower1 = rep(denseEigen.SparsePCA.err.mean - 2*denseEigen.SparsePCA.err.sd, 2),
                                         CI.upper2 = rep(denseEigen.reguPCA.err.mean + 2*denseEigen.reguPCA.err.sd, 2),
                                         CI.lower2 = rep(denseEigen.reguPCA.err.mean - 2*denseEigen.reguPCA.err.sd, 2))


qplot(SampleSize, Err.percen, data=denseEigen.SparsePCA.err.df, geom="line", col = Method,
      main="Error Percentage in Case 1-Dense Eigenvector",size=I(2))+
        geom_ribbon(aes(ymin=CI.lower1,ymax=CI.upper1), alpha=0.1,col = "darkgreen")+
        geom_ribbon(aes(ymin=CI.lower2,ymax=CI.upper2), alpha=0.1,col = "red")+
        geom_hline(yintercept=True.Error.Prop,linetype="dashed", color = "black",size=1)


                                         
        