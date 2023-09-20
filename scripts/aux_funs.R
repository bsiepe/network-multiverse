# taken from https://github.com/aweigard/GIMME_AR_simulations/blob/master/simulate_data.R

# This function generates group-level relations

# input parameters:
# p.con = proportion of paths that are contemporaneous (vs. lagged)
# nvar = number of variables in the time series
# AR = average strength of autoregressive realtions (SD=.10)
# dens = density of  (not including autoregressive relations)
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
#' @param A A matrix (for no subgroups) or list of A matrices, with slice # = # of subgroups. Contemporaneous Effects.
#' @param Phi Phi matrix (for no subgroups) or list of Phi matrices, with slice # = # of subgroups. Temporal Effects.
#' @param Psi matrix (for no subgroups) or list of Psi matrices, with slice # = # of subgroups. Prediction Errors.
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
                               split_graph = TRUE, # show temporal and contemporaneous seperated 
                               n_ind = NULL,      # specify number of individuals manually for now
                               endogeneous = FALSE,   # is there an endogeneous var? currently allows only 1 endog.
                               temp_labels = NULL, # labels for temporal network
                               cont_labels = NULL, # labels for contemporaneous network
                               # pct = FALSE,      # use percentages for plotting
                               ...){              # additional arguments to be passed to qgraph 
  # browser()
  
  # Count matrix across iterations
  # congratulations to myself for this unwieldy code
  count_mat <- Reduce('+', lapply(mv_res$l_adj_i, function(x){
    x <- Filter(function(entry) is.double(entry), x)
    as.matrix(Reduce('+', x))
  }))

  
  # n_ind <- length(mv_res$n_ind[[1]])
  

  
  # Convert to proportions
  # TODO divide by number of individuals
  prop_mat <- count_mat / (nrow(mv_res) * n_ind)
  
  # Cutoff option
  if(!is.null(cutoff)){
    prop_mat[abs(prop_mat) < cutoff] <- 0
  }
  
  # Split matrix based on lag vs. non-lagged
  if(is.null(n_lagged)){
    # round down for odd numbers
    n_lagged <-  floor(ncol(prop_mat)/2)
  }
  temp_mat <- prop_mat[, 1:n_lagged]
  cont_mat <- prop_mat[, (n_lagged + 1): ncol(prop_mat)]
  
  # Account for endogenous variable
  # insert zero row for the endogeneous variable
  if(isTRUE(endogeneous)){
    cont_mat <- rbind(cont_mat, rep(0, n_lagged + 1))
    if(is.null(temp_labels)){
      temp_labels <- colnames(cont_mat)[-length(cont_mat)]
    }
  }
  if(is.null(cont_labels)){
    cont_labels <- colnames(cont_mat)
  }
  if(is.null(temp_labels)){
    temp_labels <- colnames(cont_mat)
  }
  
  #--- Compute percentage labels
  pct_labels_temp <- paste0(round(t(temp_mat)*100,1), " %")
  
  pct_labels_cont <- paste0(round(t(cont_mat)*100,1), " %")
  
  # Plot
  if(isTRUE(split_graph)){
    # par(mfrow = c(1,2))
    # Uncomment this if layout should be done automatically
    # Temporal
    qgraph::qgraph(
      input = t(temp_mat), 
      layout       = "circle",
      lty = 1,
      edge.labels  = pct_labels_temp,
      theme = "colorblind",
      negDashed = TRUE,
      parallelEdge = TRUE,
      # arrows       = FALSE,
      labels = temp_labels,    # so that this does not show "lag" in name
      # title.cex = 2, 
      # label.cex    = 1.25,
      title = "Temporal",
      maximum = 1,
      ...)
    
    # Contemporaneous
    qgraph::qgraph(
      input = t(cont_mat), 
      layout  = "circle",
      lty = 1,
      edge.labels  = pct_labels_cont,
      theme = "colorblind",
      negDashed = TRUE,
      parallelEdge = TRUE,
      # arrows   = FALSE,
      labels = cont_labels,
      # title.cex = 2,
      # label.cex = 1.25,
      title = "Contemporaneous",
      maximum = 1, 
      ...)
    
  }
  
  
  if(isFALSE(split_graph)){
    qgraph::qgraphMixed(
      undirected = t(prop_mat[,(n_lagged+1):(ncol(prop_mat))]),
      directed = t(prop_mat[, 1: (n_lagged)]),
      layout       = "circle",
      ltyUndirected = 1,
      ltyDirected = 2,
      edge.labels  = TRUE,
      edge.color   = "blue",
      parallelEdge = TRUE,
      # arrows       = FALSE,
      labels       = 
        colnames(prop_mat)[(n_lagged+1):(ncol(prop_mat))],
      # title.cex = 2,
      # label.cex    = 1.25,
      ...)
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



# Specification Plot ------------------------------------------------------

# For lineplot
plot_outcome <- function(mv_res, 
                         var,
                         specs = NULL, # hard-coded for now
                         y_label){   
  
  var <- enquo(var)
  
  mv_res %>% 
    dplyr::mutate(variable = as.numeric(!!var)) %>% 
    dplyr::arrange(variable) %>% 
    dplyr::mutate(iteration = dplyr::row_number()) %>%
    ggplot(aes(x = .data$iteration,
               y = variable)) + 
    geom_point(size = 0.8)+
    theme_multiverse()+
    labs(x = "", 
           y = y_label)
}

# For Specification Plot
plot_specification <- function(mv_res,
                               var,
                               specs = NULL){    # hard-coded for now
  
  var <- enquo(var)
  
  mv_res %>% 
    dplyr::mutate(variable = as.numeric(!!var)) %>% 
    dplyr::arrange(variable) %>% 
    dplyr::mutate(iteration = dplyr::row_number()) %>%
    dplyr::mutate(across(c(groupcutoffs, subcutoffs,
                           rmsea.cuts, srmr.cuts,
                           cfi.cuts, nnfi.cuts,
                           n.excellent), ~ as.factor(.))) %>% 
    dplyr::mutate(groupcutoffs = fct_recode(groupcutoffs, !!!setNames(as.character(group_cuts), group_levels)),
                  subcutoffs = fct_recode(subcutoffs, !!!setNames(as.character(sub_cuts), sub_levels)),
                  rmsea.cuts = fct_recode(rmsea.cuts, !!!setNames(as.character(rmsea_cuts), rmsea_levels)),
                  srmr.cuts = fct_recode(srmr.cuts, !!!setNames(as.character(srmr_cuts), srmr_levels)),
                  cfi.cuts = fct_recode(cfi.cuts, !!!setNames(as.character(cfi_cuts), cfi_levels)),
                  nnfi.cuts = fct_recode(nnfi.cuts, !!!setNames(as.character(nnfi_cuts), nnfi_levels)),
                  n.excellent = fct_recode(n.excellent, !!!setNames(as.character(n_excels), n_excels_levels))) %>% 
    tidyr::pivot_longer(cols = c(groupcutoffs,
                                 subcutoffs,
                                 rmsea.cuts,
                                 srmr.cuts,
                                 cfi.cuts,
                                 nnfi.cuts,
                                 n.excellent),
                        values_to = "value", names_to = "specification") %>%
    dplyr::mutate(specification = dplyr::case_match(specification,
                                                    "groupcutoffs" ~ "Group",
                                                    "subcutoffs" ~ "Subgroup",
                                                    "rmsea.cuts" ~ "RMSEA",
                                                    "srmr.cuts" ~ "SRMR",
                                                    "cfi.cuts" ~ "CFI",
                                                    "nnfi.cuts" ~ "NNFI",
                                                    "n.excellent" ~ "N° excellent")) %>% 
    dplyr::mutate(specification = as.factor(specification)) %>% 
    dplyr::mutate(specification = forcats::fct_relevel(specification, 
                                                       "Group", 
                                                       "Subgroup",
                                                       "RMSEA", 
                                                       "SRMR",
                                                       "CFI",
                                                       "NNFI")) %>% 
    dplyr::mutate(value = forcats::fct_relevel(value, 
                                               "liberal",
                                               "medium-liberal",
                                               "medium",
                                               "medium-strict",
                                               "strict")) %>% 
    ggplot(aes(x = .data$iteration,
               y = 1,
               color = .data$value)) + 
    geom_point(shape = 124, size = 15
               #pch='.'   #for faster plotting
    )+
    scale_y_continuous(limits = c(0.99, 1.01), expand = c(0,0))+
    theme_multiverse()+
    scale_color_manual(values = palette_full)+
    facet_wrap(specification~., 
               ncol = 1, 
               strip.position = "left")+
    labs(y = "",
         x = "Iteration",
         color = "Specification")+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          legend.text = element_text(size = rel(1.3)),
          legend.title = element_text(size = rel(1.3)),
          legend.spacing.x = unit(0.2, "cm"),
          legend.text.align = 0,
          strip.text = element_text(size = rel(1.4)),
          axis.text.x = element_text(size = rel(1.2)),
          axis.title.x = element_text(size = rel(1.3)))
}


# -------------------------------------------------------------------------


# MVGIMME STUFF -----------------------------------------------------------

# -------------------------------------------------------------------------


#' Compare Results from Multiverse Analysis
#'
#' This function compares the results obtained from a multiverse analysis with a reference model. It performs three levels of comparison: group-level, subgroup-level, and individual-level comparison. The results are combined into a tibble where each row represents one specification, and individual results can be stored in a dataframe inside the tibble.
#'
#' @param l_res A list of results obtained from the multiverse analysis. Each element of the list should be a data object containing the results.
#' @param ref_model The reference model to compare the results with.
#' @return A tibble containing the comparison results at different levels: group-level, subgroup-level, and individual-level comparison along with condition information for each specification.
#' @export
multiverse.compare <- function(l_res,
                               ref_model
){
  # browser()
  
  #---- Compare Group-Level
  comp_group <- multiverse.compare.group(
    l_res = l_res, ref_model = ref_model
  )
  
  #---- Compare Subgroup-Level
  comp_subgroup <- multiverse.compare.subgroup(
    l_res = l_res, ref_model = ref_model
  )  
  
  #---- Compare Individual Solutions
  comp_ind <- multiverse.compare.individual(
    l_res = l_res, ref_model = ref_model
  )   
  
  #--- Output
  # Idea: tibble where each row is one specification
  # individual results can be stored in dataframe inside tibble
  comp_res <- dplyr::bind_cols(comp_group, comp_subgroup, comp_ind)
  
  
  # save condition for each specification 
  conds <- do.call(rbind, lapply(l_res, `[[`, "conds"))
  comp_res <- cbind(comp_res, conds)
  
  return(comp_res)
}


multiverse.compare.group <- function(l_res,
                                     ref_model){
  
  #--- Reference model info
  n_var <- ref_model$n_vars_total
  
  # indices for temporal and contemporaneous
  temp_ind <- 1:(n_var/2)
  cont_ind <- ((n_var/2)+1):n_var
  
  #--- Adjacency Matrix
  ## Find group effects adjacency matrix
  # find n converged for reference model
  n_ind_ref <- sum(unlist(lapply(ref_model$path_est_mats, is.double)))
  
  # Find converged for multiverse
  n_ind <- lapply(l_res, function(x){
    sum(unlist(lapply(x$path_est_mats, is.double)))
  })
  
  ref_groupedge <- ifelse(ref_model$path_counts == n_ind_ref, 1, 0)
  # ignore autoregressive coefs
  diag(ref_groupedge[, temp_ind]) <- rep(0, n_var/2)
  
  # Find group effects adjacency matrix differences
  l_ref_diff <- lapply(l_res, function(x){
    n_ind_mv <- sum(unlist(lapply(x$path_est_mats, is.double)))
    tmp_groupedge <- ifelse(x$path_counts == n_ind_mv, 1, 0)
    diag(tmp_groupedge[, temp_ind]) <- rep(0, n_var/2)
    diff_groupedge <- ref_groupedge - tmp_groupedge
    return(diff_groupedge)
  }
  )
  
  # Count occurrence of each group effect
  ## list of adjacency matrices
  l_adjacency <- lapply(l_res, function(x){
    n_ind_mv <- sum(unlist(lapply(x$path_est_mats, is.double)))
    tmp_groupedge <- ifelse(x$path_counts == n_ind_mv, 1, 0)
    diag(tmp_groupedge[, temp_ind]) <- rep(0, n_var/2)
    return(tmp_groupedge)
  }
  )
  
  #--- Heterogeneity
  ## divide no. of group edges by no. of total edges
  l_heterogeneity <- list()
  for(i in 1:length(l_adjacency)){
    # calculate number of estimated edges, group + individual
    tmp_mat <- l_res[[i]]$path_counts
    n_ind_mv <- sum(unlist(lapply(l_res[[i]]$path_est_mats, is.double)))
    diag(tmp_mat[, temp_ind]) <- rep(0, n_var/2)
    l_heterogeneity[[i]] <- sum(tmp_mat[tmp_mat == n_ind_mv]) / sum(tmp_mat)
  }
  
  #--- Output
  l_out <- tibble(
    n_ind = n_ind,
    ref_diff_g = l_ref_diff, 
    adjacency_g = l_adjacency,
    heterogeneity_g = l_heterogeneity
  ) %>% 
    dplyr::mutate(ref_diff_g = as.list(ref_diff_g)) %>% 
    tidyr::unnest(n_ind) %>% 
    tidyr::unnest(heterogeneity_g)
  
  return(l_out)
}


multiverse.compare.subgroup <- function(l_res, 
                                        ref_model){
  
  #--- Reference model info
  n_ind <- length(ref_model$data)
  ref_sim_matrix <- ref_model$sim_matrix
  
  #--- Number of subgroups
  l_n_sub <- lapply(l_res, function(x){
    return(length(unique(x$fit$sub_membership)))
  }
  )
  
  
  #--- Size of subgroups
  l_size_sub <- lapply(l_res, function(x){
    return(table(x$fit$sub_membership))
  })
  
  #--- Subgroup similarity matrix & modularity
  ## Use similarity matrix of each specification
  l_pert <- lapply(l_res, function(x){
    l_calc <- list()
    # can use VI and ARI from perturbR package here
    l_calc$vi <- vi.dist(ref_sim_matrix, x$sim_matrix)
    l_calc$ari <- arandi(ref_sim_matrix, x$sim_matrix)
    l_calc$modularity <- as.numeric(x$fit$modularity[1])   
    return(l_calc)
  })
  
  #--- Subgroup edges
  ## Loop over subgroups
  l_out <- tibble(
    n_sub_g = l_n_sub,
    size_sub_s = l_size_sub,
    l_pert = l_pert
  ) %>% 
    tidyr::unnest_wider(l_pert)
  
  return(l_out)
}



multiverse.compare.individual <- function(l_res,
                                          ref_model){
  
  # browser()
  
  #--- Reference model 
  n_ind <- length(ref_model$path_est_mats)
  # n_var <- ref_model$n_lagged + ref_model$n_endog
  # the former version did not include TimeOfDay
  n_var <- ref_model$n_vars_total
  
  # indices for temporal and contemporaneous
  temp_ind <- 1:(n_var/2)
  cont_ind <- ((n_var/2)+1):n_var
  
  ## Estimates
  ref_path_est_mats <- ref_model$path_est_mats
  
  
  #--- Adjacency matrix
  # ignore autoregressive effects
  ref_adj_mats <- lapply(ref_path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      tmp <- ifelse(x != 0, 1, 0)
      diag(tmp[, temp_ind]) <- rep(0, n_var/2)
      return(tmp)
    }
  })
  
  #--- Density
  ref_dens_temp <- lapply(ref_path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      abs_sum(x[, temp_ind])
    }
  })
  
  ref_dens_cont <- lapply(ref_path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      abs_sum(x[, cont_ind])
    }
  })
  
  #--- Fit indices
  fit_ind_names <- c("chisq", "df", "npar", "pvalue", "rmsea", "srmr",
                     "nnfi", "cfi", "bic", "aic", "logl")
  ref_fit_ind <- ref_model$fit[names(ref_model$fit) %in% fit_ind_names]
  
  #--- Centrality
  ref_outstrength <- lapply(ref_path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      colSums(abs(x))
    }
  })
  
  ########################
  #--- Nonconverg. Checks
  ########################
  # TODO
  # can only add this after I have seen some actual nonconvergence
  
  
  ########################
  #--- Plausibility Checks
  ########################
  
  # Check diagonal of Psi (contemporaneous) matrix
  # Values should be >= 0 & <= 1
  # inspired by https://github.com/aweigard/GIMME_AR_simulations/blob/master/analyze_recovery_Balanced.R
  # 1 means implausible value somewhere in psi matrix
  # l_implausible <- lapply(l_res, function(x){
  #   unlist(lapply(x$path_est_mats, function(y){
  #     if(!is.double(y)){
  #       NA
  #     }
  #     else{
  #       sum(ifelse(any(diag(y[,temp_ind]) < 0 | diag(y[,temp_ind]) > 1) , 1, 0))
  #     }
  #     
  #   }))
  # })
  # sum_implausible <- sapply(l_implausible, sum)   # this might throw an error if parts are NA
  
  
  ########################
  #--- Compare Edges
  ########################
  #--- Nondirectional recovery
  # Only for contemporaneous effects
  # ref_nondir_adj_mats <- lapply(ref_path_est_mats, function(x){
  #   tmp <- ifelse(x[,cont_ind] != 0, 1, 0)
  #   tmp <- nondirect_adjacency(tmp)
  #   return(tmp)
  # })
  # 
  # l_diff_nondir_adj <- lapply(l_res, function(x){
  #   tmp_nondir_adj_mats <- lapply(x$path_est_mats, function(y){
  #     tmp <- ifelse(y[,cont_ind] != 0, 1, 0)
  #     tmp <- nondirect_adjacency(tmp)
  #     return(tmp)
  #   })
  #   l_nondir_adj <- Map('-', ref_nondir_adj_mats, tmp_nondir_adj_mats)
  #   lapply(l_nondir_adj, function(y){as.matrix(y)})
  #   
  # })
  
  #--- Directional recovery
  l_diff_adj <- lapply(l_res, function(x){
    ## adjacency matrices
    # ignore AR coefs
    tmp_adj_mats <- lapply(x$path_est_mats, function(y){
      if(!is.double(y)){
        NA
      }
      else{
        tmp <- ifelse(y != 0, 1, 0)
        diag(tmp[, temp_ind]) <- rep(0, n_var/2)
        return(tmp)
      }

    })
    l_adj <- Map('-', ref_adj_mats, tmp_adj_mats)
    lapply(l_adj, function(y){as.matrix(y)})
  })
  
  # Return adjacency matrix for each individual
  l_adj <- lapply(l_res, function(x){
    tmp_adj_mats <- lapply(x$path_est_mats, function(y){
      if(!is.double(y)){
        NA
      }
      else{
        tmp <- ifelse(y != 0, 1, 0)
        diag(tmp[, temp_ind]) <- rep(0, n_var/2)
        return(tmp)
      }

    })
    return(tmp_adj_mats)
  })
  
  #--- Difference/bias path estimates
  l_diff_ests <- lapply(l_res, function(x){
    l_est <- Map('-', ref_path_est_mats, x$path_est_mats)
    l_est <- lapply(l_est, function(y){
      as.matrix(y)})
    # filter entries that are NA
    l_est <- Filter(function(entry) is.double(entry), l_est)
    simplify2array(l_est)
  })
  
  
  #--- Density
  l_diff_dens_temp <- lapply(l_res, function(x){
    tmp_denstemp <- lapply(x$path_est_mats, function(y){
      if(!is.double(y)){
        NA
      }
      else{
        abs_sum(y[, temp_ind])        
      }
    })
    unlist(Map('-', ref_dens_temp, tmp_denstemp))
  })
  
  l_diff_dens_cont <- lapply(l_res, function(x){
    tmp_denscont <- lapply(x$path_est_mats, function(y){
      if(!is.double(y)){
        NA
      }
      else{
        abs_sum(y[, cont_ind])
      }
    })
    unlist(Map('-', ref_dens_cont, tmp_denscont))
  })
  
  #--- Compare Fits
  l_diff_fit <- lapply(l_res, function(x){
    return(ref_fit_ind - x$fit[names(x$fit) %in% fit_ind_names])
  })
  
  
  #--- Compare centrality
  # compute outstrength
  l_cent <- lapply(l_res, function(x){
    tmp_outstrength <- lapply(x$path_est_mats, function(y){
      if(!is.double(y)){
        NA
      }
      else{
        colSums(abs(y))
      }
      
    })
  })
  
  # compute difference
  l_diff_cent <- lapply(l_res, function(x){
    tmp_outstrength <- lapply(x$path_est_mats, function(y){
      if(!is.double(y)){
        NA
      }
      else{
        colSums(abs(y))
      }
    })
    Map('-', ref_outstrength, tmp_outstrength)
  })
  
  ## binary: most central node the same?
  # split by temporal and contemporaneous
  central_node_identical <- list()
  
  for(i in seq_along(l_cent)){
    central_node_identical[[i]] <- list()
    for(j in seq_along(ref_outstrength)){
      max_temp_ref <- which.max(ref_outstrength[[j]][temp_ind])
      max_temp_mv <- which.max(l_cent[[i]][[j]][temp_ind])
      max_cont_ref <- which.max(ref_outstrength[[j]][cont_ind])
      max_cont_mv <- which.max(l_cent[[i]][[j]][cont_ind])
      central_node_identical[[i]][[j]] <- list()
      central_node_identical[[i]][[j]]$temp_identical <- if(!is.null(max_temp_ref) & !is.null(max_temp_mv)){
        names(max_temp_ref) == names(max_temp_mv)}
      else{
        NA
      }
      central_node_identical[[i]][[j]]$cont_identical <- if(!is.null(max_cont_ref) & !is.null(max_cont_mv)){
        names(max_cont_ref) == names(max_cont_mv)}
      else{
        NA
      }
    }
  }
  
  
  ########################
  #--- Aggregate
  ########################
  # First aggregate by taking mean for each edge across individuals
  # then summarize this again
  # split by averaging over all nonzero differences, or all differences
  
  #--- Nondirected adjacency
  # mean_diff_nondir_adj <- lapply(l_diff_nondir_adj, function(x){
  #   l_tmp <- list()
  #   l_tmp$diff_nondir_adj_sum_mat_i <- apply(simplify2array(x), 1:2, abs_sum)
  #   l_tmp$diff_nondir_adj_sum_sum_i <- sum(l_tmp$nondir_adj_sum_mat)
  #   l_tmp$diff_nondir_adj_sum_mean_i <- mean(l_tmp$nondir_adj_sum_mat)
  #   return(l_tmp)
  # })
  
  #--- Adjacency matrix
  mean_diff_adj <- lapply(l_diff_adj, function(x){
    l_tmp <- list()
    # Filter entries that are NA (and therefore saved as 1x1 integer)
    x <- Filter(function(entry) is.double(entry), x)
      l_tmp$diff_adj_sum_mat_i <- apply(simplify2array(x), 1:2, abs_sum)
      l_tmp$diff_adj_sum_sum_i <- sum(l_tmp$diff_adj_sum_mat_i, na.rm = TRUE)
      l_tmp$diff_adj_sum_mean_i <- mean(l_tmp$diff_adj_sum_mat_i, na.rm = TRUE)

    return(l_tmp)
  })
  
  
  #--- Mean differences of edges
  mean_diff_ests <- lapply(l_diff_ests, function(x){
    l_tmp <- list()
    # mean_mat <- apply(simplify2array(x), 1:2, mean)
    l_tmp$mean_nonzero_diff_edge_i <- abs_mean(x[which(x != 0, arr.ind = TRUE)])
    l_tmp$med_nonzero_diff_edge_i <- abs_med(x[which(x != 0, arr.ind = TRUE)])
    l_tmp$mean_diff_edge_i <- abs_mean(x)
    l_tmp$med_diff_edge_i <- abs_med(x)
    return(l_tmp)
  })
  
  #--- Densities
  mean_abs_diff_dens_temp <- sapply(l_diff_dens_temp, abs_mean)
  mean_abs_diff_dens_cont <- sapply(l_diff_dens_cont, abs_mean)
  mean_diff_dens_temp <- sapply(l_diff_dens_temp, mean, na.rm = TRUE)
  mean_diff_dens_cont <- sapply(l_diff_dens_cont, mean, na.rm = TRUE)
  
  
  
  #--- Fits 
  mean_diff_fits <- lapply(l_diff_fit, function(x){
    l_tmp <- list()
    l_tmp$mean_diff_fit_i <- apply(x, 2, abs_mean)
    l_tmp$med_diff_fit_i <- apply(x, 2, abs_med)
    return(l_tmp)
  })
  
  
  #--- Centrality
  # difference across all centrality values
  mean_diff_cent <- lapply(l_diff_cent, function(x){
    l_tmp <- list()
    # filter out NA
    x <- Filter(function(entry) is.double(entry), x)
    diff_cent <- rowSums(simplify2array(x), na.rm = TRUE)
    l_tmp$mean_diff_cent_i <- abs_mean(diff_cent)
    l_tmp$med_diff_cent_i <- abs_med(diff_cent)
    return(l_tmp)
  })
  
  # central node identical
  sum_temp_central_identical <- lapply(central_node_identical, function(x) 
    lapply(x, function(y) 
      sum(unlist(y$temp_identical), na.rm = TRUE)))
  sum_cont_central_identical <- lapply(central_node_identical, function(x) 
    lapply(x, function(y) 
      sum(unlist(y$cont_identical), na.rm = TRUE)))
  
  sum_temp_central_identical <- sapply(sum_temp_central_identical, function(x) sum(unlist(x), na.rm = TRUE))
  sum_cont_central_identical <- sapply(sum_cont_central_identical, function(x) sum(unlist(x), na.rm = TRUE))
  
  
  l_out <- tibble(
    # l_implausible_i = l_implausible,
    # sum_implausible_i = sum_implausible,
    # l_diff_nondir_adj_i = l_diff_nondir_adj,
    l_adj_i = l_adj, 
    l_diff_adj_i = l_diff_adj,
    l_diff_ests_i = l_diff_ests,
    l_diff_fit_i = l_diff_fit,
    l_diff_cent_i = l_diff_cent,
    central_node_identical_i = central_node_identical,
    # mean_diff_nondir_adj_i = mean_diff_nondir_adj,
    mean_diff_adj_i = mean_diff_adj,
    mean_diff_ests_i = mean_diff_ests,
    mean_diff_cent_i = mean_diff_cent,
    mean_abs_diff_dens_temp_i = mean_abs_diff_dens_temp,
    mean_abs_diff_dens_cont_i = mean_abs_diff_dens_cont,
    mean_diff_dens_temp_i = mean_diff_dens_temp,
    mean_diff_dens_cont_i = mean_diff_dens_cont,
    sum_temp_central_identical_i = sum_temp_central_identical,
    sum_cont_central_identical_i = sum_cont_central_identical,
    mean_diff_fit_i = mean_diff_fits
  ) %>% 
    tidyr::unnest_wider(c(
                          # mean_diff_nondir_adj_i,
                          # central_node_identical_i,
                          mean_diff_adj_i,
                          mean_diff_ests_i, 
                          mean_diff_cent_i,
                          mean_diff_fit_i))
  
  return(l_out)
}








# From perturbr internal functions ----------------------------------------
# https://github.com/cran/perturbR/blob/master/R/vi.dist.R
vi.dist <-
  function(cl1,cl2,parts=FALSE, base=2){ 
    if(length(cl1) != length(cl2)) stop("cl1 and cl2 must have same length")
    
    # entropy 
    ent <- function(cl){
      n <- length(cl)
      p <- table(cl)/n
      -sum(p*log(p, base=base))
    }
    # mutual information
    mi <- function(cl1,cl2){
      p12 <- table(cl1,cl2)/length(cl1)
      p1p2 <- outer(table(cl1)/length(cl1),table(cl2)/length(cl2))
      sum(p12[p12>0]*log(p12[p12>0]/p1p2[p12>0], base=base))
    }
    
    if(!parts) return(ent(cl1) + ent(cl2) -2*mi(cl1,cl2))
    ent1 <- ent(cl1)
    ent2 <- ent(cl2)
    mi12 <- mi(cl1,cl2)
    c("vi"=ent1+ent2-2*mi12, "H(1|2)" =ent1-mi12, "H(2|1)"=ent2 -mi12)
  }

# https://github.com/cran/perturbR/blob/master/R/arandi.R
arandi <-
  function(cl1,cl2, adjust=TRUE){
    if(length(cl1)!=length(cl2)) stop("cl1 and cl2 must have same length")
    tab.1 <- table(cl1)
    tab.2 <- table(cl2)
    tab.12 <- table(cl1,cl2)
    if(adjust){
      correc <- sum(choose(tab.1,2))*sum(choose(tab.2,2))/choose(length(cl2),2)
      return((sum(choose(tab.12,2))-correc)/(0.5*sum(choose(tab.1,2))+0.5*sum(choose(tab.2,2))-correc) )}
    else{ 
      1+(sum(tab.12^2)-0.5*sum(tab.1^2)-0.5*sum(tab.2^2))/choose(length(cl2),2)
    }
  }



# Matrix list summary -----------------------------------------------------
matrix_summary <- function(l_matrix) {
  # Calculate mean, median, sd, min, and max for each matrix element across the list
  arr_list <- simplify2array(l_matrix)
  
  a_mean <- apply(arr_list, 1:2, mean)
  a_median <- apply(arr_list, 1:2, stats::median)
  a_min <- apply(arr_list, 1:2, min)
  a_max <- apply(arr_list, 1:2, max)
  a_sd <- apply(arr_list, 1:2, sd)
  
  # Combine the results into a data frame
  summary_df <- as.data.frame(summary_list)
  
  return(summary_df)
}


# Small helpers -----------------------------------------------------------

abs_mean <- function(x){
  mean(abs(x), na.rm = TRUE)
}
abs_med <- function(x){
  stats::median(abs(x), na.rm = TRUE)
}

abs_sum <- function(x){
  sum(abs(x), na.rm = TRUE)
}


nondirect_adjacency <- function(adj_mat) {
  # Number of ariables
  n_adj_vars <- nrow(adj_mat)
  
  # Initialize symmetrical matrix with 0s
  sym_matrix <- matrix(0, nrow = n_adj_vars, ncol = n_adj_vars)
  
  # Iterate through each cell of the original matrix
  for (i in 1:n_adj_vars) {
    for (j in 1:n_adj_vars) {
      # If there is any effect (1) in either direction, update the symmetrical matrix
      if (adj_mat[i, j] == 1 | adj_mat[j, i] == 1) {
        sym_matrix[i, j] <- 1
        sym_matrix[j, i] <- 1
      }
    }
  }
  
  return(sym_matrix)
}








