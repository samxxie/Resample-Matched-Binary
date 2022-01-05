
## R code for simulating covid social data (type I error comparison, change of statistic)

net_xxie_dir <- "C:\\Users\\xxie\\Dropbox (EinsteinMed)\\ODrivePartial\\Xianhong Xie\\Mimi\\DOM\\Covid_Social\\Data";
#net_xxie_dir <- "\\\\tsclient\\C\\Users\\xxie\\Dropbox (EinsteinMed)\\ODrivePartial\\Xianhong Xie\\Mimi\\DOM\\Covid_Social\\Data";
setwd(net_xxie_dir)

fixefp1 <- c(-2.37, 0.32, -0.97, 0.30, -0.18, -0.75, -0.75, -0.81, -0.61, -1.72) # intercept followed by 9 question effects (q1 is ref)
#ranefp1 <- c(0.09, 1.22) # estimated CS parameters from the glimmix model
#ranefp1 <- c(0.60, 0.70) # correlation damped G matrix
ranefp1 <- c(0.30, 1.00) # correlation damped G matrix
#ranefp1 <- c(0.10, 1.20) # correlation damped G matrix

n <- 300
gmat <- ranefp1[1]*diag(10) + ranefp1[2]*rep(1,10)%*%t(rep(1,10))
cholg <- chol(gmat)
X <- cbind(1, kronecker(matrix(1,n,1),kronecker(matrix(1,2,1),diag(10))[,-1]))
Z <- kronecker(matrix(1,2*n,1), diag(10))

iseed1 <- 100
iseed2 <- 200

nsimul <- 500
#nsimul <- 10000
test_res_ovl <- rep(NA, nsimul)
test_res_ind <- matrix(NA, nsimul, 10)
colnames(test_res_ind) <- paste("sdh", 1:10, sep="")

test_res_ovlp <- test_res_ovl
test_res_indp <- test_res_ind

#samp_method  <- "Boot1"
##samp_method <- "Boot2"
#samp_method <- "Permute1"
#samp_method <- "Permute2"
samp_method <- "BootPerm"
nboot <- 10000


## define the functions

get.isimul.dat <- function(isimul, iseed, X, Z, fixefp, cholg) {
  
  set.seed(iseed)
  
  n <- nrow(X) / 20
  sdat <- data.frame(matrix(NA, 2*n, 12))
  names(sdat) <- c("id", "pd", paste("sdh",1:10,sep=""))
  
  if (isimul > 1) {
    for (j in 1:(isimul-1)) {
      for (i in 1:n) {
        tmp1 <- sum(rnorm(10))
        tmp2 <- sum(sapply(1:20, function(i) rbinom(1,1,0.5))) 
      }
    }
  }
  
  for (i in 1:n) {
    indexi <- ((i-1)*20+1):(i*20)
    Xi <- X[indexi,]
    Zi <- Z[indexi,]
    etveci  <- Xi%*%fixefp + Zi%*%as.vector(t(cholg)%*%rnorm(10))
    sdhveci <- sapply(1/(1+exp(-etveci)), function(p) rbinom(1,1,p))
    sdat[((i-1)*2+1):(i*2),] <- cbind(i, 1:2, matrix(sdhveci,2,10,byrow=T))
  }
  
  sdat
  
}


#apply(sdat1[,-c(1:2)], 2, mean)
#cor(sdat1[,-c(1:2)])

#apply(sdat1[sdat1$pd==1,-c(1:2)], 2, mean)
#apply(sdat1[sdat1$pd==2,-c(1:2)], 2, mean)


rank_chg_test2 <- function(sdat, iseed=200, nboot=10000, samp_method="Boot1") {
  
  meanSum1b  <- t(apply(sdat[,-(1:2)], 2, function(x) tapply(x,sdat$pd,mean)))
  meanRank1b <- apply(-meanSum1b, 2, rank)
  
  meanRank2b <- cbind(meanRank1b, diff=apply(meanRank1b,1,function(x) x[1]-x[2]))
  #summary(meanRank2b[,3])
  #table(meanRank2b[,3])
  
  samplepd1b <- sdat$pd
  
  set.seed(iseed)
  size1b  <- nrow(sdat)
  #stat0b  <- max(abs(meanRank2b[,3]))
  stat0b  <- sum(abs(meanRank2b[,3]))
  stat0jb <- meanRank2b[,3]
  statv1b <- rep(NA, nboot)
  statm1b <- matrix(NA, nboot, 10)
  colnames(statm1b) <- rownames(meanRank2b)
  #samp_method <- "Boot1"
  #samp_method <- "Boot2"
  #samp_method <- "Permute1"
  #samp_method <- "Permute2"
  #samp_method <- "BootPerm"
  for (i in 1:nboot) {
    
    if (samp_method=="Boot1") {
      meanSum  <- t(apply(sdat[sample(size1b,repl=T),-(1:2)], 2, function(x) tapply(x,samplepd1b,mean)))
    } else if (samp_method=="Boot2") {
      meanSum  <- t(apply(sdat[as.vector(sapply(0:(size1b/2-1),function(i) 2*i+sample(2,repl=T))),-(1:2)], 
                  2, function(x) tapply(x,samplepd1b,mean)))
    } else if (samp_method=="Permute1") { # ** Note: permute at subject level does not change the rank stats **
      meanSum  <- t(apply(sdat[sample(size1b,repl=F),-(1:2)], 2, function(x) tapply(x,samplepd1b,mean)))
    } else if (samp_method=="Permute2") {
      meanSum  <- t(apply(sdat[as.vector(sapply(0:(size1b/2-1),function(i) 2*i+sample(2,repl=F))),-(1:2)], 
                  2, function(x) tapply(x,samplepd1b,mean)))
    } else if (samp_method=="BootPerm") {
      meanSum  <- t(apply(sdat[as.vector(sapply(0:(size1b/2-1),function(i) 2*sample(0:(size1b/2-1),1,repl=T)+sample(2,repl=F))),-(1:2)], 
                  2, function(x) tapply(x,samplepd1b,mean)))
    }
    
    meanRank <- apply(-meanSum, 2, rank) 
    
    rankDiff <- apply(meanRank,1,function(x) x[1]-x[2]) 
    
    #statv1b[i]  <- max(abs(rankDiff))
    statv1b[i]  <- sum(abs(rankDiff))
    statm1b[i,] <- rankDiff
    
    if ((i%%500)==0) { 
      cat(paste("  Iteration",i,"done\n",sep=" ")) 
    }
  }
  
  pval1b <- mean(stat0b<=statv1b)
  pval1jb <- apply(rbind(stat0jb,statm1b), 2, function(x) { 
    if (x[1] > 0) { 
      p1 <- mean(x[-1]>=x[1]) 
    } else if (x[1] < 0) {
      p1 <- mean(x[-1]<=x[1]) 
    } else {
      p1 <- min(mean(x[-1]>=x[1]), mean(x[-1]<=x[1]))
    }
    
    min(2*p1,1)
  })
  pval1mb <- cbind("Rank before"=meanRank1b[,1], "Rank after"=meanRank1b[,2], "Rank change"=stat0jb, "P-value"=pval1jb)
  
  list(pval.ovl = pval1b, pval.ind=pval1mb)
  
}


## run the simulations

#for (isimul in 1:1) {
for (isimul in 1:nsimul) {
  cat("Working on simulation", isimul, "\n")
  
  sdat1 <- get.isimul.dat(isimul, iseed1, X, Z, fixefp1, cholg)
  #rank_test_res <- rank_chg_test(sdat1, iseed=iseed2, nboot=nboot, samp_method=samp_method)
  rank_test_res <- rank_chg_test2(sdat1, iseed=iseed2, nboot=nboot, samp_method=samp_method)
  
  test_res_ovl[isimul]  <- (rank_test_res$pval.ovl < 0.05)*1
  test_res_ind[isimul,] <- (rank_test_res$pval.ind[,4] < 0.05)*1
  
  test_res_ovlp[isimul]  <- rank_test_res$pval.ovl
  test_res_indp[isimul,] <- rank_test_res$pval.ind[,4]
  
  cat("\n")
}

list(n=n, nsimul=nsimul, samp_method=samp_method)
mean(test_res_ovl, na.rm=T)
apply(test_res_ind, 2, mean, na.rm=T)
write.csv(t(c(mean(test_res_ovl,na.rm=T),matrix(apply(test_res_ind,2,mean,na.rm=T),1,10))*100), file="..\\Results\\tmp.csv", quote=F, row.names=F)

mean(test_res_ovlp < 0.05, na.rm=T)
apply(test_res_indp, 2, function(x) mean(x < 0.05,na.rm=T))

siglev <- sort(unique(test_res_ovlp))
t1err  <- sapply(sort(unique(test_res_ovlp)), function(z) mean(test_res_ovlp < z, na.rm=T))

#plot(siglev, t1err, type="l", xlim=c(0,1), ylim=c(0,1))
#lines(c(0,1), c(0,1), lty=2)
##abline(v=0.05, lty=3)


#save.image(file="..\\Results\\type1_error_method_boot1_result_sz300_ns500s.RData")
##save.image(file="..\\Results\\type1_error_method_boot2_result_sz300_ns500s.RData")
#save.image(file="..\\Results\\type1_error_method_permute1_result_sz300_ns500s.RData")
#save.image(file="..\\Results\\type1_error_method_permute2_result_sz300_ns500s.RData")
save.image(file="..\\Results\\type1_error_method_bootperm_result_sz300_ns500s.RData")


###############################################################################
