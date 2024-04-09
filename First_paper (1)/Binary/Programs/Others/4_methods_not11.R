library(doParallel)
library(clue)
library(matrixStats)
library(klaR)
library(ludic)
#library(doBy)

# Fellegi-Sunter with 3 categorical comparison
FS3 <- function(datA, datB, K){
  nA = nrow(datA)
  nB = nrow(datB)
  comp_mat <- compare_binary1(datA, datB, K=K)
  
  ## Using fixed starting point
  fit = EM_1(comp_mat, datA, datB, K, tol=1e-6, maxits = 500)
  g = fit$g
 
  threshold1 = 0.5
  TPR1 = sum(comp_mat[which(g>threshold1), K+1])/sum(comp_mat[, K+1])
  PPV1 = sum(comp_mat[which(g>threshold1), K+1])/length(which(g>threshold1))
  
  #threshold2 = 0.9
  #TPR2 = sum(comp_mat[which(g>threshold2), K+1])/sum(comp_mat[, K+1])
  #PPV2 = sum(comp_mat[which(g>threshold2), K+1])/length(which(g>threshold2))
  
  return(c(TPR_1, PPV1))
}


###### Boris
boris <- function(datA, datB, K){
  nB = nrow(datB)
  nA = nrow(datA)
  boris = recordLink(datB[,1:K], datA[,1:K], eps_plus =0.01, eps_minus = 0.01,use_diff = FALSE)
  
  threshold = 0.5
  upper = which(boris > threshold)
  indexB = upper%%nB
  indexB[which(indexB==0)] = nB
  indexA = ceiling(upper/nB)
  idB = datB[indexB, K+1]
  idA = datA[indexA, K+1]
  
  TPR_boris_threshold = sum(idA==idB)/nB
  PPV_boris_threshold = sum(idA==idB)/length(idB)
  
  return(c(TPR_boris_11, TPR_boris_threshold, PPV_boris_threshold))
}

FS <- function(datA, datB, K){
  nA = nrow(datA)
  nB =  nrow(datB)
  
  comp_mat <- compare_binary(datA, datB, K)
  
  fit <- EM_binary_pfix(comp_mat, datA, datB, K,tol = 1e-6, maxits = 300)
  g = fit$g
  g[g<0]=0
  Gmat = matrix(g, nrow = nB, byrow = TRUE)
  opti_binary <- solve_LSAP(Gmat, maximum = TRUE)
  predict_binary = cbind(seq_along(opti_binary), opti_binary)
  
  # Percentage of correct link
  TPR_FS11 = sum(predict_binary[,2] == datB[,K+1])/nB
  
  threshold = 0.5
  upper = which(Gmat >= threshold)
  indexB = upper%%nB
  indexB[which(indexB==0)] = nB
  indexA = ceiling(upper/nB)
  idB = datB[indexB, K+1]
  idA = datA[indexA, K+1]
  
  TPR_FS_threshold = sum(idA==idB)/nB
  PPV_FS_threshold = sum(idA==idB)/length(idB)
  return(c(TPR_FS11,TPR_FS_threshold,PPV_FS_threshold))
}



### Fellegi-Sunter with 4 categorical comparison
FS_2_fixed <- function(datA, datB, K){
  nB = nrow(datB)
  comp_mat2 <- compare_binary2(datA, datB, K=K)
  
  ## Using fixed starting point
  fit_2 = EM_2(comp_mat2, datA, datB, K, tol=1e-5, maxits = 500)
  g = fit_2$g
  g[g<0]=0
  Gmat = matrix(g, nrow = nB, byrow = TRUE)
  opti_cate <- solve_LSAP(Gmat, maximum = TRUE)
  predict_cate = cbind(seq_along(opti_cate), opti_cate)
  # Percentage of correct link
  TPR_2= sum(predict_cate[,2] == datB[,K+1])/nB
  #Using threshold
  threshold = 0.5
  upper = which(Gmat >= threshold)
  indexB = upper%%nB
  indexB[which(indexB==0)] = nB
  indexA = ceiling(upper/nB)
  idB = datB[indexB, K+1]
  idA = datA[indexA, K+1]
  
  TPR_2_threshold = sum(idA==idB)/nB
  PPV_2_threshold = sum(idA==idB)/length(idB)
  
  
  ## Using fixed parameters  
  fit_2 = EM_2_fixed(comp_mat2, datA, datB, K, tol=1e-5, maxits = 500)
  g = fit_2$g
  g[g<0]=0
  Gmat = matrix(g, nrow = nB, byrow = TRUE)
  opti_cate <- solve_LSAP(Gmat, maximum = TRUE)
  predict_cate = cbind(seq_along(opti_cate), opti_cate)
  # Percentage of correct link
  TPR_2_fixed= sum(predict_cate[,2] == datB[,K+1])/nB
  
  threshold = 0.5
  upper = which(Gmat >= threshold)
  indexB = upper%%nB
  indexB[which(indexB==0)] = nB
  indexA = ceiling(upper/nB)
  idB = datB[indexB, K+1]
  idA = datA[indexA, K+1]
  
  TPR_2_threshold_fixed = sum(idA==idB)/nB
  PPV_2_threshold_fixed = sum(idA==idB)/length(idB)
  
  ## Using true estimates
  fit_true = EM_true_2(comp_mat2, datA, datB, K)
  g = fit_true$g
  g[g<0]=0
  Gmat = matrix(g, nrow = nB, byrow = TRUE)
  opti_cate <- solve_LSAP(Gmat, maximum = TRUE)
  predict_cate = cbind(seq_along(opti_cate), opti_cate)
  # Percentage of correct link
  TPR_true= sum(predict_cate[,2] == datB[,K+1])/nB
  
  
  
  return(c(TPR_2, TPR_2_fixed,TPR_true, TPR_2_threshold, PPV_2_threshold, TPR_2_threshold_fixed, PPV_2_threshold_fixed))
}

