# add simulateVARtest

simulateVARtest <- function(A         = NULL, 
                            Phi       = NULL, 
                            Psi       = NULL, 
                            subAssign = NULL, 
                            N         = NULL, 
                            ASign     = "random",  
                            PhiSign   = "random",  
                            tp       = NULL,
                            indA      = 0.01, 
                            indPhi    = 0.01,
                            indPsi    = 0.00) {
  
  # browser()
  
  AList   <- list()
  PhiList <- list()
  PsiList <- list()
  dataList  <- list()
  subgroups <- unique(subAssign)
  
  ###### Checks and Rules #####
  
  if(!is.list(Phi)){
    PhiList[[1]] <- Phi
  }
  
  # BS: Added this: If Phi is a List, PhiList is empty now. This throws an error
  # 
  if(is.list(Phi)){
    PhiList <- Phi
  }
  
  if(is.null(Psi)){
    Psi <- matrix(0, nrow = dim(PhiList[[1]])[1], ncol = dim(PhiList[[1]])[2])
    diag(Psi) <- 1
  }
  
  if(is.null(A)){
    A <- matrix(0, nrow = dim(PhiList[[1]])[1], ncol = dim(PhiList[[1]])[2])
  }
  
  if(is.list(A)) {
    AList   <- A
    PhiList <- Phi
    PsiList <- Psi
    if(is.null(subAssign)){
      writeLines("NOTICE: Multiple patterns provided with no subgroup assignments",
                 "By default subgroup assignments generated with equal prob. for each ind.")
      nSubs <- length(A)
      subAssign <- matrix(,N, 1)
      for (ind in 1:N)
        subAssign[ind] <- sample((1:nSubs),1)
    }
    if(length(is.list(A)) != length(is.list(Phi))) 
      stop(paste0("ERROR: Different numbers of matrices provided for A and Phi.",
                  " Please ensure there is one A and Phi matrix for each subgroup."))
    if(length(is.list(A)) != length(is.list(Psi))) 
      stop(paste0("ERROR: Different numbers of matrices provided for A and Psi.",
                  " Please ensure there is one A and Psi matrix for each subgroup."))
  } else {
    writeLines("NOTICE: One A matrix provided. No subgroups generated.")
    subAssign <- matrix(1, N, 1)
    AList[[1]]   <- A
    PsiList[[1]] <- Psi
  }
  
  if(is.null(N))
    stop(paste0("ERROR: Please provide N for the number of individuals to generate"))
  if(is.null(tp))
    stop(paste0("ERROR: Please provide tp for the number of time points per person"))
  
  ##### Data Generation ####
  
  vars <- length(AList[[1]][1,])
  
  Ind.nonstation <- c()
  
  # storing individual parameter matrices
  indList <- list()
  
  for (ind in 1:N) {
    go <- 1
    counter <- 0                                   
    while (go > 0 & counter < 100){
      ATemp <- AList[[subAssign[ind]]]
      PhiTemp <- PhiList[[subAssign[ind]]]
      PsiTemp <- PsiList[[subAssign[ind]]]
      
      AMean <- mean(ATemp)
      PhiMean <- mean(PhiTemp)
      PsiMean <- mean(PsiTemp)
      
      # BS: change such that indA applies to individual elements,
      # not all elements of that person 
      
      
      
      if(indA>0){
        n.col.a <- ncol(ATemp)
        n.row.a <- nrow(ATemp)
        
        # Loop over individual elements
        for(i in 1:n.col.a){
          for(j in 1:n.row.a){
            if(ATemp[i,j] == 0){
              ATemp[i,j] <- stats::rbinom(1, 1, indA)
            }
            if(ATemp[i,j] == 1){
              if(ASign == "random"){
                random.sign <- sample(c(0,1), size=1)
                if(random.sign==0) {ATemp[i,j] <- -stats::rnorm(1, AMean, 0.3)}
                if(random.sign==1) {ATemp[i,j] <- stats::rnorm(1, AMean, 0.3)}
              }
              if(ASign == "neg")
                ATemp[i,j] <- -stats::rnorm(1, AMean, 0.3)
              if(ASign == "pos")
                ATemp[i,j] <- stats::rnorm(1, AMean, 0.3)
            }  
          }
        }
      }
      diag(ATemp) <- 0
      
      
      if(indPhi>0){
        n.col.phi <- ncol(ATemp)
        n.row.phi <- nrow(ATemp)
        
        # Loop over individual elements
        for(i in 1:n.col.phi){
          for(j in 1:n.row.phi){
            if(PhiTemp[i,j] == 0){
              PhiTemp[i,j] <- stats::rbinom(1, 1, indPhi)
            }
            if(PhiTemp[i,j] == 1){
              if(PhiSign == "random"){
                random.sign <- sample(c(0,1), size=1)
                if(random.sign==0) {PhiTemp[i,j] <- -stats::rnorm(1, PhiMean, 0.3)}
                if(random.sign==1) {PhiTemp[i,j] <- stats::rnorm(1, PhiMean, 0.3)}
              }
              if(PhiSign == "neg")
                PhiTemp[i,j] <- -stats::rnorm(1, PhiMean, 0.3)
              if(PhiSign == "pos")
                PhiTemp[i,j] <- stats::rnorm(1, PhiMean, 0.3)
            }  
          }
        }
      }
      
      if(indPsi>0){
        PsiTemp[which(PsiTemp == 0)] <- stats::rbinom(1,1,indPsi)
        PsiTemp[which(PsiTemp == 1)] <- stats::rnorm(1, PsiMean, 0.3)
      }
      
      # BS: Return individual parameter lists
      indList[[ind]] <- list()
      # indList[[ind]]$AInd <- ATemp
      # indList[[ind]]$PhiInd <- PhiTemp
      indList[[ind]]$PsiInd <- PsiTemp
      indList[[ind]]$Paths <- cbind(PhiTemp, ATemp)
      # combine into path counts matrix
      indList[[ind]]$Adj <- as.matrix(ifelse(cbind(PhiTemp, ATemp) != 0, 1, 0))
      
      negA <- solve(diag(vars)-ATemp) 
      
      time  <- matrix(0, nrow = vars, ncol = tp+400) # 400 burn-in observations
      
      time1 <- matrix(0, nrow = vars, ncol = tp+400)
      
      
      # BS: Gates et al. 2017, Eq. 3
      noise <- negA %*% t(MASS::mvrnorm(n = (tp+400),rep(0,vars),PsiTemp, empirical = TRUE))
      # BS: this throws an non-conformable error if PsiTemp is not properly defined
      
      time[,1] <- noise[,1]
      
      time1[,1] <- negA %*% PhiTemp %*% time[,1] + noise[,1]
      
      time[,2]  <- time1[,1]
      
      # BS: Gates et al. 2017, Eq. 3
      for (t in 2:(tp+400)){
        time1[,t]  <- negA %*% PhiTemp %*% time[,(t)] + noise[,t]
        if (t<(tp+400))
          time[,(t+1)] <- time1[,t]
      }
      go <- 0
      for (c in 1:length(time[,1])){
        adf_result <- suppressWarnings(tryCatch({aTSA::adf.test(time[c, ], out = FALSE)},
                                                error = function(e) NA))       # BS: changed this
        if(adf_result$type3[1,3]>0.05 || 
           is.na(adf_result$type3[1,3])
           ||
           sum(abs(time[c, 400:(tp+400)])) > 10000)    # BS: add check for large values of timeseries
          go <- go + 1
        counter <- sum(counter, 1)
      }
    }
    if(counter == 100){
      Ind.nonstation <- append(Ind.nonstation, ind)
      writeLines(paste0('WARNING: No Stationary Time Series Data Generated for Individual:', ind))
    } else {
      dataList[[ind]] <- t(time[,401:(400+tp)]) 
      names(dataList)[ind]<- paste0('ind', ind)
    }
  }
  
  
  return(list(A=AList,
              Phi = PhiList, 
              Psi = PsiList, 
              dataList = dataList, 
              subAssign = subAssign,
              indList = indList))
  
}







# Summarize helpers -------------------------------------------------------

#--- Compute Precision
precision <- function(true, est){
  # true paths relative to all paths
  true_pos <- sum(true & est)
  all_pos <- sum(est)
  
  if (all_pos == 0) {
    return(0)
  } else {
    precision <- true_pos / all_pos
    return(precision)
  }
}

#--- Compute Recall
recall <- function(true, est) {
  true_pos <- sum(true & est)
  all_true <- sum(true)
  
  if (all_true == 0) {
    return(0)
  } else {
    recall <- true_pos / all_true
    return(recall)
  }
}




#--- Nondirected Adjacency matrix
nondirect_adjacency <- function(adj_mat) {
  # Number of ariables
  n_adj_vars <- nrow(adj_mat)
  
  # Initialize symmetrical matrix with 0s
  sym_matrix <- matrix(0, nrow = n_adj_vars, ncol = n_adj_vars)
  
  # Iterate through each cell of the original matrix
  for (i in 1:n_adj_vars) {
    for (j in 1:n_adj_vars) {
      # If there is any effect (1) in either direction, update the symmetrical matrix
      if (n_adj_vars[i, j] == 1 || n_adj_vars[j, i] == 1) {
        sym_matrix[i, j] <- 1
        sym_matrix[j, i] <- 1
      }
    }
  }
  
  return(sym_matrix)
}


#--- Some absolute summary stats
abs_mean <- function(x){
  mean(abs(x), na.rm = TRUE)
}

abs_med <- function(x){
  stats::median(abs(x), na.rm = TRUE)
}

abs_sum <- function(x){
  sum(abs(x), na.rm = TRUE)
}
