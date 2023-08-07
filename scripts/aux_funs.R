# taken from https://github.com/aweigard/GIMME_AR_simulations/blob/master/simulate_data.R

# This function generates group-level relations

# input parameters:
# p.con = proportion of paths that are contemporaneous (vs. lagged)
# nvar = number of variables in the time series
# AR = average strength of autoregressive realtions (SD=.10)
# dens = density of network (not including autoregressive relations)
# p.group = proportion of relations at the group (vs. individual) level
# con.b = average strength of contemporaneous relations (SD=.10)
# lag.b = average strength of cross-lagged relations (SD=.10)


mat.generate.asw <- function(p.con, nvar,AR,dens,p.group,con.b,lag.b){
  repeat{
    A   <- matrix(0, ncol = nvar, nrow = nvar) # contemporaneous
    Phi <- matrix(0, ncol = nvar, nrow = nvar) # lagged (AR in diag)
    pos.all <- 2*(nrow(A)*ncol(A)) - 2*nrow(A) # all available (non-ar) spots
    cnt.all <- dens*p.group*pos.all # number of group paths given density and proportion of group paths
    indices <- which(Phi == 0, arr.ind = TRUE)
    indices <- indices[which(indices[,1] != indices[,2]), ]# kick out diag
    row.col <- sample(1:nrow(indices), cnt.all, replace = F)# sample n=cnt.all paths 
    # group paths (in case there are an odd number, randomly chose whether con or lag gets one extra)
    n.p.1<-row.col[1:round(p.con*length(row.col))]
    n.p.2<-row.col[(length(n.p.1)+1):length(row.col)]
    n.p<-list(n.p.1,n.p.2)
    rand<-sample(1:2,size=2,replace=FALSE)
    grp.con <- n.p[[ rand[1] ]]
    grp.lag <- n.p[[ rand[2] ]]
    
    Phi[indices[grp.lag,]] <- lag.b
    A[indices[grp.con,]]   <- con.b
    diag(Phi)              <- AR # Insert AR terms here!
    
    # disallow bidirectional contemporaneous paths and A matrices with max eigenvalues >1
    A.test<-1*(A!=0)
    A.test<-(A.test+t(A.test)) 
    if ( (max(A.test)!=2) & (max(abs(eigen(A, only.values = FALSE)$values))<1) ) break
    
  }
  
  all <- cbind(Phi, A)
  ind.pres <- which(all != 0, arr.ind = T)
  ind.pres <- ind.pres[which(ind.pres[,1] != ind.pres[,2]), ]
  
  level <- "grp"
  
  all.lvl           <- matrix(NA, ncol = ncol(all), nrow = nrow(all))
  all.lvl[ind.pres] <- level
  diag(all.lvl)     <- "grp"
  
  all_sub1 <- all
  all_lvl1 <- all.lvl
  
  res <- list(sub1 = all_sub1, 
              lvl1 = all_lvl1)
  return(res)
}



# This function generates group-level relations adds individual-level
# relations to the matrix, adds noise, and simulates the time series

# input parameters:
# mat = matrix of group-level relations
# lvl = matric of input relation levels (e.g., "group")
# p.con = proportion of paths that are contemporaneous (vs. lagged)
# p.group = proportion of relations at the group (vs. individual) level
# con.b = average strength of contemporaneous relations (SD=.10)
# lag.b = average strength of cross-lagged relations (SD=.10)

ts.generate.asw <- function (mat, lvl, t,dens,p.group,con.b,lag.b,p.con) {
  repeat {
    
    repeat{
      v <- ncol(mat)/2 #calc nvars from matrix
      Phi <- mat[, 1:v] # pull out Phi
      A   <- mat[, (v+1):(v*2)] # pull out A
      A_ind        <- matrix(0, ncol=v, nrow=v)
      ## finds indices of zero elements in matrix
      indices.A      <- which(A==0,arr.ind=T)
      indices.Phi      <- which(Phi==0,arr.ind=T)
      ## removes indices of diagonal elements from consideration
      indices.A      <- indices.A[which(indices.A[,1]!=indices.A[,2]),]
      indices.Phi      <- indices.Phi[which(indices.Phi [,1]!=indices.Phi[,2]),]
      
      ## determines the number of individual paths to add
      # (in case there are an odd number, randomly chose whether con or lag gets one extra)
      pos.all<-(v*v-v)*2
      cnt.all <- dens*p.group*pos.all  
      rand<-sample(c(round(cnt.all/2),(round(cnt.all)-round(cnt.all/2))),size=2,replace=FALSE)
      ## randomly selects row/col combinations for individual-level paths
      row.col.A      <- sample(1:nrow(indices.A), rand[1], replace = F) 
      row.col.Phi      <- sample(1:nrow(indices.Phi), rand[2], replace = F) 
      
      # Betas for lagged and contemporaneous paths
      Phi[indices.Phi[row.col.Phi,]] <- lag.b
      A[indices.A[row.col.A,]]     <- con.b
      
      # add noise to A betas, SD =.1
      noise.inds      <- which(A != 0, arr.ind = TRUE) 
      A[noise.inds]   <- A[noise.inds] + rnorm(n = nrow(noise.inds), mean = 0, sd = .1)
      
      # add noise to Phi betas, excluding ar terms, SD = .1
      noise.inds      <- which(Phi != 0, arr.ind = TRUE) 
      noise.inds      <- noise.inds[which(noise.inds[,1] != noise.inds[,2]), ]# kick out diag
      Phi[noise.inds] <- Phi[noise.inds] + rnorm(n = nrow(noise.inds), mean = 0, sd = .1)
      
      # add noise to AR terms, SD = .1 (if you want)
      noise.inds      <- which(Phi !=Inf, arr.ind = TRUE)    
      noise.inds      <- noise.inds[which(noise.inds[,1] == noise.inds[,2]), ]# include diag
      Phi[noise.inds] <- Phi[noise.inds] + rnorm(n = nrow(noise.inds), mean = 0, sd = .1)
      
      #disallow bidirectional contemporaneous paths and and A matrices with max eigenvalues >1
      A.test<-1*(A!=0)
      A.test<-(A.test+t(A.test)) 
      if( (max(A.test)!=2) & (max(abs(eigen(A, only.values = FALSE)$values))<1) ) break
    }
    
    
    st <- (t+50) #Alex added, now robust to any t 
    noise <- matrix(rnorm(v*st,0,1),v) #
    I     <- diag(v) # identity matrix
    time  <- matrix(0,nrow=v, ncol=(st+1))
    time1 <- matrix(0,nrow=v, ncol=st)
    
    # simulate data points for each time step
    for (i in 1:st){
      time1[,i]  <- solve(I-A)%*%(Phi%*%time[,i] + noise[,i])
      time[,i+1] <- time1[,i]
    }               
    time1  <- time1[,(51:(50+t))] # Fixed this, was initially 50:, causing ts to be one too long
    series <- t(time1)
    paths  <- cbind(Phi, A)
    if (abs(max(series, na.rm = TRUE)) < 20 & abs(min(series, na.rm = TRUE)) > .01 
        & abs(min(series, na.rm = TRUE)) < 20) break
  }
  
  lvl[is.na(lvl) & paths != 0] <- "ind"
  
  list   <- list("series"  = series,
                 "paths"   = paths,
                 "levels"  = lvl)
  return(list)
}






# New function that combines both -----------------------------------------
sim.gimme <- function(p.con, 
                      nvar,
                      AR,
                      dens,
                      p.group,
                      con.b,
                      lag.b,
                      mat, 
                      lvl,         # 
                      t,           # length of time series
                      n.ind        # number of individuals
                      ){
  
  # Group-level
  group.mat <- mat.generate.asw(
                   p.con = p.con,
                   nvar = nvar,
                   AR = AR,
                   dens = dens,
                   p.group = p.group,
                   con.b = con.b,
                   lag.b = lag.b) 
  
  
  # Individual level
  # loop over number of individuals 
  data.list <- list()
  for(i in 1:n.ind){
    ind.data <- ts.generate.asw(
      mat = group.mat$sub1,
      lvl = group.mat$lvl1,
      t = t,
      p.group = p.group,
      con.b = con.b,
      lag.b = lag.b,
      p.con = p.con,
      dens = dens
      
    )
    ind.data$series <- round(ind.data$series, digits = 5)
  
    data.list[[i]] <- ind.data
    }

  return(data.list)
  
  
} 










# New simulation function from gimme --------------------------------------

#' @name simulateVAR
#' @title Simulate data from Vector AutoRegression (VAR) models.
#' @description This function simulates data. It allows for structural VAR and VAR data generating models. 
#' @usage
#' simulateVAR(A   = NULL, 
#'             Phi       = NULL, 
#'             Psi       = NULL, 
#'             subAssign = NULL, 
#'             N         = NULL, 
#'             ASign     = "random",  
#'             PhiSign   = "random",  
#'             Obs       = NULL,
#'             indA      = 0.01, 
#'             indPhi    = 0.01,
#'             indPsi    = 0.00)
#' @param A A matrix (for no subgroups) or list of A matrices, with slice # = # of subgroups. 
#' @param Phi Phi matrix (for no subgroups) or list of Phi matrices, with slice # = # of subgroups.  
#' @param Psi matrix (for no subgroups) or list of Psi matrices, with slice # = # of subgroups. 
#' @param subAssign Optional vector of length N that indicates which subgroup each individual is in. 
#' @param N Number of indvidiuals.
#' @param Obs Number of observations (T) per individual. Burn in of 400 is used to generate then discarded.  
#' @param indA Sparsity of individual-level A paths. 0 indicates no individual-level. Use decimals. Default is 
#' 0.01, meaning that each path that is not in the group-level A matrix has a 0.01 chance of being added. 
#' @param indPhi Sparsity of individual-level Phi paths. 0 indicates no individual-level. Use decimals. Default is 
#' 0.01, meaning that each path that is not in the group-level Phi matrix has a 0.01 chance of being added.
#' @param indPsi Sparsity of individual-level Psi paths. 0 indicates no individual-level. Use decimals. Default is 
#' 0, meaning that each path that is not in the group-level Psi matrix has a 0 chance of being added at the ind. level.
#' Individual- level paths added at this rate per individual. 
#' @param ASign Defaults to "random" for ind level paths, with 50 percent chance of positive and 50 percent negative, other option is either "neg" or "pos" which provides all negative or all positive relations, respectively.
#' @param PhiSign Defaults to "random" for ind level paths, with 50 percent chance of positive and 50 percent negative, other option is either "neg" or "pos" which provides all negative or all positive relations, respectively. 
#' @author KM Gates, Ai Ye, Ethan McCormick, & Zachary Fisher 
#' @export simulateVAR

simulateVARtest <- function(A         = NULL, 
                        Phi       = NULL, 
                        Psi       = NULL, 
                        subAssign = NULL, 
                        N         = NULL, 
                        ASign     = "random",  
                        PhiSign   = "random",  
                        Obs       = NULL,
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
  if(is.null(Obs))
    stop(paste0("ERROR: Please provide Obs for the number of time points per person"))
  
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
      
      time  <- matrix(0, nrow = vars, ncol = Obs+400) # 400 burn-in observations
      
      time1 <- matrix(0, nrow = vars, ncol = Obs+400)
      
      noise <- negA %*% t(MASS::mvrnorm(n = (Obs+400),rep(0,vars),PsiTemp, empirical = TRUE))
      # BS: this throws an non-conformable error if PsiTemp is not properly defined
      
      time[,1] <- noise[,1]
      
      time1[,1] <- negA %*% PhiTemp %*% time[,1] + noise[,1]
      
      time[,2]  <- time1[,1]
      
      for (t in 2:(Obs+400)){
        time1[,t]  <- negA %*% PhiTemp %*% time[,(t)] + noise[,t]
        if (t<(Obs+400))
          time[,(t+1)] <- time1[,t]
      }
      go <- 0
      for (c in 1:length(time[,1])){
        adf_result <- suppressWarnings(tryCatch({aTSA::adf.test(time[c, ], out = FALSE)},
                               error = function(e) NA))       # BS: changed this
        if(adf_result$type3[1,3]>0.05 || 
           is.na(adf_result$type3[1,3])
           ||
           sum(abs(time[c, 400:(Obs+400)])) > 10000)    # BS: add check for large values of timeseries
          go <- go + 1
        counter <- sum(counter, 1)
      }
    }
    if(counter == 100){
      Ind.nonstation <- append(Ind.nonstation, ind)
      writeLines(paste0('WARNING: No Stationary Time Series Data Generated for Individual:', ind))
    } else {
      dataList[[ind]] <- t(time[,401:(400+Obs)]) 
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

























# -------------------------------------------------------------------------
# VISUALIZATION -----------------------------------------------------------
# -------------------------------------------------------------------------

# Multiverse Network Plot -------------------------------------------------
# This function sums all adjacency matrices of all individuals
# across all specifications, and then creates a plot where
# thickness of graphs corresponds to number of inclusions (irrespective of sign of edge)
# see here: https://github.com/GatesLab/gimme/blob/cb0cf2f6b1cf5db5b16330966ccd8920cef15c66/gimme/R/summaryPathsCounts.R#L139

multiverse.network <- function(mv_res, 
                               n_lagged = NULL, # number of lagged variables, assumed to be half the columns if not specified
                               cutoff = NULL,    # only include effects with certain proportion of occurrence? 
                               split_graph = TRUE){  # show temporal and contemporaneous seperated    
  
  # Count matrix across iterations
  count_mat <- as.matrix(Reduce('+', mv_res$adj_sum_mat_i))
  
  

  
  # Convert to proportions
  prop_mat <- count_mat / nrow(mv_res)
  
  # Cutoff option
  if(!is.null(cutoff)){
    prop_mat[abs(prop_mat) < cutoff] <- 0
  }
  
  # Split matrix based on lag vs. non-lagged
  if(is.null(n_lagged)){
    n_lagged <-  ncol(prop_mat)/2
  }
  temp_mat <- prop_mat[, 1:n_lagged]
  cont_mat <- prop_mat[, (n_lagged + 1): ncol(prop_mat)]
  
  # Plot
  if(isTRUE(split_graph)){
    par(mfrow = c(1,2))
    
    # Temporal
    qgraph::qgraph(
      input = temp_mat, 
      layout       = "circle",
      lty = 1,
      edge.labels  = TRUE,
      theme = "colorblind",
      negDashed = TRUE,
      parallelEdge = TRUE,
      fade         = FALSE,
      # arrows       = FALSE,
      labels = colnames(cont_mat),    # so that this does not show "lag" in name
      label.cex    = 2,
      title = "Temporal")
    
    # Contemporaneous
    qgraph::qgraph(
      input = cont_mat, 
      layout  = "circle",
      lty = 1,
      edge.labels  = TRUE,
      theme = "colorblind",
      negDashed = TRUE,
      parallelEdge = TRUE,
      fade = FALSE,
      # arrows   = FALSE,
      labels = colnames(cont_mat),
      label.cex = 2,
      title = "Contemporaneous")
    
  }
  
  
  if(isFALSE(split_graph)){
    qgraph::qgraphMixed(
      undirected = prop_mat[,(n_lagged+1):(ncol(prop_mat))],
      directed = prop_mat[, 1: (n_lagged)],
      layout       = "circle",
      ltyUndirected = 1,
      ltyDirected = 2,
      edge.labels  = TRUE,
      edge.color   = "blue",
      parallelEdge = TRUE,
      fade         = FALSE,
      # arrows       = FALSE,
      labels       = 
        colnames(prop_mat)[(n_lagged+1):(ncol(prop_mat))],
      label.cex    = 2)
  }

  
}



# Detrending --------------------------------------------------------------
fn_detrend <- function(x,                    
                       vars,                 
                       time_var = "time",   
                       sig_only = FALSE){    
  for (v in 1:length(vars)){
    # Regress on time
    lm_form <- as.formula(paste0(vars[v], "~", time_var))
    lm_res <- summary(lm(lm_form, data = x))
    # detrend with residuals
    # [,4] accesses p-values
    # [2] p-value of beta of tp
    if(sig_only){
      if(lm_res$coefficients[,4][2] < 0.05){
        x[!is.na(x[vars[v]]),vars[v]] <- residuals(lm_res)
      }
    }
    if(isFALSE(sig_only)){
      print(paste0("Detrend variable ", vars[v]))
      x[!is.na(x[vars[v]]),vars[v]] <- residuals(lm_res)
    }
    
  }
  return(x)
}



# ggplot2 theme -----------------------------------------------------------
theme_multiverse <- function() {
  # add google font
  sysfonts::font_add_google("News Cycle", "news")
  # use showtext
  showtext::showtext_auto()
  # theme
  ggplot2::theme_minimal(base_family = "news") +
    ggplot2::theme(
      # remove minor grid
      panel.grid.minor = ggplot2::element_blank(),
      # Title and Axis Texts
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2), hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1.1), hjust = 0.5),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, b = 10)),
      
      # Faceting
      strip.text = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.1), hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      # Grid
      panel.grid = ggplot2::element_line(colour = "#F3F4F5"),
      # Legend
      legend.title = ggplot2::element_text(face = "plain"),
      legend.position = "top",
      legend.justification = 1,
      # Panel/Facets
      panel.spacing.y = ggplot2::unit(1.5, "lines")
    )
}




