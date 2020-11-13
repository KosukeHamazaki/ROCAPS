require(rrBLUP)
require(BGLR)
source(paste0(dir.script, "1.1.0_Rice_ROCAPS_MTM2.R"))


MTM_prediction <- function(Y, Y.test.no = NULL, XF = NULL, x = NULL, Z = NULL, K = NULL, n.lines = NULL, NDCG.k = 10,
                           resCov = list(type = "UN", df0 = ncol(as.matrix(Y)), S0 = diag(ncol(as.matrix(Y)))), env.num = 1,
                           nIter = 110, burnIn = 10, thin = 2, saveAt = "", tolD = 1e-05, Parallel = TRUE, n.core = NULL,
                           CV = TRUE, n.fold = 10, seed = NULL, silent = FALSE, time = TRUE, verbose = FALSE){
  st0 <- Sys.time()
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  #### Some calculation for kinship matrix ####
  if(is.null(K)){
    if(is.null(x)){
      stop("x or K must be non-null!!!")
    }
    K <- A.mat(x)
  }
  
  if(is.null(n.lines)){
    n.lines <- nrow(K)
  }
  
  if(!is.null(Z)){
    K <- tcrossprod(Z %*% K, Z) 
  }
  
  #### Some modification of phenotype data ####
  Y <- as.matrix(Y)
  num.pheno <- ncol(Y)
  
  
  if(CV){
    #### n-fold cros-svalidation ####
    id <- sample(1:n.lines %% n.fold)
    id[id == 0] <- n.fold
    
    ids <- rep(id, env.num)
    
    
    y.obs <- Y[, ncol(Y)]
    y.pred <- rep(NA, nrow(Y))
    
    
    
    #### The start of cross-validation ####
    if(Parallel){
      if(is.null(n.core)){
        .cores <- parallel::detectCores()
      }else{
        if(n.core <= 0) {
          .cores <- parallel::detectCores() + n.core
        } else{
          .cores <- n.core
        }
      }
      doParallel::registerDoParallel(.cores)
      fms <- foreach(i = 1:n.fold, .export = ls(envir = parent.frame()),
                     .packages = "BGLR") %dopar% {
                       
                       #### Prepare the test data-set ####
                       Y.test <- Y
                       Y.test[ids == i, ncol(Y)] <- NA
                       
                       
                       #### The start of prediction by MTM ####
                       saveAt.fold <- paste0(saveAt, "/fold_", i, "/")
                       dir.create(saveAt.fold)
                       if(num.pheno >= 2){
                         fm <- MTM2(
                           Y = Y.test,
                           XF = XF,
                           K = list(
                             list(
                               K = K,
                               COV = list(
                                 type = 'UN',
                                 df0 = num.pheno,
                                 S0 = diag(num.pheno)
                               )
                             )
                           ),
                           resCov = resCov,
                           nIter = nIter,
                           burnIn = burnIn,
                           thin = thin,
                           saveAt = saveAt.fold,
                           tolD = tolD,
                           verbose = verbose
                         )
                       }else{
                         ETA <- list(K = list(K = K, model = "RKHS"))
                         fm <- BGLR(
                           y = Y.test,
                           ETA = ETA,
                           nIter = nIter,
                           burnIn = burnIn,
                           thin = thin,
                           saveAt = saveAt.fold,
                           verbose = verbose
                         ) 
                       }
                       return(fm)
                     }
      
      for(i in 1:n.fold){
        if(num.pheno >= 2){
          y.pred[ids == i] <- (fms[[i]])$YHat[ids == i, ncol(Y)]
        } else {
          y.pred[ids == i] <- (fms[[i]])$yHat[ids == i]
        }
      }
    } else {
      fms <- NULL
      
      #### The start of cross-validation ####
      for(i in 1:n.fold){
        #### Display the iteration no. of the cross-validation ####
        if(!silent){
          if(i != 1){
            print(paste0("n.fold ", i, "/", n.fold, " : It took ", end - st, " ", attr(end - st, "units"), " for this validation."))
          }else{
            print(paste0("n.fold ", i, "/", n.fold, " : This is the start of cross-validation.")) 
          }
        }
        
        st <- Sys.time()
        
        
        #### Prepare the test data-set ####
        Y.test <- Y
        Y.test[ids == i, ncol(Y)] <- NA
        
        
        #### The start of prediction by MTM ####
        saveAt.fold <- paste0(saveAt, "/fold_", i, "/")
        dir.create(saveAt.fold)
        if(num.pheno >= 2){
          fm <- MTM2(
            Y = Y.test,
            XF = XF,
            K = list(
              list(
                K = K,
                COV = list(
                  type = 'UN',
                  df0 = num.pheno,
                  S0 = diag(num.pheno)
                )
              )
            ),
            resCov = resCov,
            nIter = nIter,
            burnIn = burnIn,
            thin = thin,
            saveAt = saveAt.fold,
            tolD = tolD,
            verbose = verbose
          )
          y.pred[ids == i] <- fm$YHat[ids == i, ncol(Y)]
        }else{
          ETA <- list(K = list(K = K, model = "RKHS"))
          fm <- BGLR(
            y = Y.test,
            ETA = ETA,
            nIter = nIter,
            burnIn = burnIn,
            thin = thin,
            saveAt = saveAt.fold,
            verbose = verbose
          ) 
          
          y.pred[ids == i] <- fm$yHat[ids == i]
        }
        fms <- c(fms, list(fm))
        end <- Sys.time()
      }
    }
  }else{
    #### When you don't perform cross-validation ####
    if(is.null(Y.test.no)){
      stop("You should assign the test data-set no. to Y.test.no when you don't perform cross-validation!")
    }
    
    y.obs <- Y[Y.test.no, ncol(Y)]
    y.pred <- rep(NA, length(Y.test.no))
    
    Y.test <- Y
    Y.test[Y.test.no, ncol(Y)] <- NA
    
    #### The start of prediction by MTM ####
    if(num.pheno >= 2){
      fm <- MTM2(
        Y = Y.test,
        XF = XF,
        K = list(
          list(
            K = K,
            COV = list(
              type = 'UN',
              df0 = num.pheno,
              S0 = diag(num.pheno)
            )
          )
        ),
        resCov = resCov,
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        saveAt = saveAt,
        tolD = tolD,
        verbose = verbose
      )
      y.pred <- fm$YHat[Y.test.no, ncol(Y)]
    }else{
      ETA <- list(K = list(K = K, model = "RKHS"))
      fm <- BGLR(
        y = Y.test,
        ETA = ETA,
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        saveAt = saveAt,
        verbose = verbose
      ) 
      y.pred <- fm$yHat[Y.test.no]
    }
    fms <- fm
  }
  
  
  if(env.num == 1){
    #### Calculate the coefficient of determination R2 ####
    R2 <- cor(y.obs, y.pred, use = "complete.obs") ^ 2
    if((n.fold == nrow(Y)) & (R2 >= 0.999)){
      R2 <- 0
    }
    
    
    #### Calculate MSE ####
    MSE <- mean((y.pred - y.obs) ^ 2, na.rm = TRUE)
    
    #### Calculate NDCG ####
    y.obs.ord <- order(y.obs, decreasing = T)
    y.pred.ord <- order(y.pred, decreasing = T)
    
    DCG.obs <- sum((y.obs[y.obs.ord])[1:NDCG.k] / log2(1 + 1:NDCG.k))
    DCG.pred <- sum((y.obs[y.pred.ord])[1:NDCG.k] / log2(1 + 1:NDCG.k))
    
    NDCG <- DCG.pred / DCG.obs
  }else{
    #### Calculate the coefficient of determination R2 ####
    R2s <- MSEs <- NDCGs <- rep(NA, env.num)
    for(i in 1:env.num){
      y.obs.now <- y.obs[((i - 1) * n.lines + 1):(i * n.lines)]
      y.pred.now <- y.pred[((i - 1) * n.lines + 1):(i * n.lines)]
      R2.now <- cor(y.obs.now, y.pred.now, use = "complete.obs") ^ 2
      
      if((n.fold == nrow(Y)) & (R2.now >= 0.999)){
        R2.now <- 0
      }
      
      #### Calculate MSE ####
      MSE.now <- mean((y.pred.now - y.obs.now) ^ 2, na.rm = TRUE)
      
      #### Calculate NDCG ####
      y.obs.ord.now <- order(y.obs.now, decreasing = T)
      y.pred.ord.now <- order(y.pred.now, decreasing = T)
      
      DCG.obs.now <- sum((y.obs.now[y.obs.ord.now])[1:NDCG.k] / log2(1 + 1:NDCG.k))
      DCG.pred.now <- sum((y.obs.now[y.pred.ord.now])[1:NDCG.k] / log2(1 + 1:NDCG.k))
      
      NDCG.now <- DCG.pred.now / DCG.obs.now
      
      
      
      R2s[i] <- R2.now
      MSEs[i] <- MSE.now 
      NDCGs[i] <- NDCG.now
    }
    
    R2 <- mean(R2s)
    MSE <- mean(MSEs)
    NDCG <- mean(NDCGs)
  }
  
  
  
  end0 <- Sys.time()
  if(time){
    print(paste0("It took ", round(end0- st0, 3), " ", attr(end0- st0, "units"), " for this prediction."))
  }
  
  return(list(fm = fms, y.obs = y.obs, y.pred = y.pred, R2 = R2, MSE = MSE, NDCG = NDCG))
}




MTM_prediction_iter <- function(Y, Y.test.no = NULL, XF = NULL, x = NULL, Z = NULL, K = NULL, n.iter = 10, fitting.scale = "R2",
                                resCov = list(type = "UN", df0 = 0, S0 = diag(0, ncol(as.matrix(Y)))), n.lines = NULL,
                                nIter = 110, burnIn = 10, thin = 2, tolD = 1e-05, Parallel = TRUE, n.core = NULL, NDCG.k = 10, env.num = 1,
                                CV = TRUE, n.fold = 10, seeds = NULL, silent = TRUE, time = TRUE, saveAt = NULL, verbose = FALSE){
  if(!is.null(seeds)){
    if(length(seeds) != n.iter){
      stop("Please set so that the length of seeds argument is same as n.iter argument!")
    }
  }
  R2s <- MSEs <- NDCGs <- rep(NA, n.iter)
  for(i in 1:n.iter){
    if(!is.null(saveAt)){
      dir.create(paste0(saveAt, "/BGLR_files_", i))
      saveAt2 <- paste0(saveAt, "/BGLR_files_", i, "/")
    }else{
      dir.create(paste0("BGLR_files_", i))
      saveAt2 <- paste0("BGLR_files_", i, "/")
    }
    res.now <- MTM_prediction(Y = Y, Y.test.no = Y.test.no, XF = XF, x = x, Z = Z, K = K, n.lines = n.lines,
                              resCov = resCov, nIter = nIter, burnIn = burnIn, thin = thin, Parallel = Parallel, n.core = n.core,
                              saveAt = saveAt2, tolD = 1e-05, NDCG.k = NDCG.k, CV = CV, env.num = env.num,
                              n.fold = n.fold, seed = seeds[i], silent = silent, time = time, verbose = verbose)
    
    
    if(!is.null(saveAt)){
      results_name <- paste0(saveAt, "/pred_results_", i, ".RData")
      save(res.now, file = results_name)
    }
    
    R2s[i] <- res.now$R2
    MSEs[i] <- res.now$MSE
    NDCGs[i] <- res.now$NDCG
  }
  
  if(fitting.scale == "R2"){
    return(R2s)
  }
  
  if(fitting.scale == "MSE"){
    return(MSEs)
  }
  
  if(fitting.scale == "NDCG"){
    return(NDCGs)
  }
  
}

