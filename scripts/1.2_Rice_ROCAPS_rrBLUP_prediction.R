require(rrBLUP)

rrBLUP_prediction <- function(Y, Y.test.no = NULL, XF = NULL, x = NULL, Z = NULL, K = NULL, n.lines = NULL, pred.method = "OCR2",
                              method = "REML", bounds = c(1e-09, 1e+09), SE = TRUE, return.Hinv = TRUE, env.num = 1,
                              NDCG.k = 10, CV = TRUE, n.fold = 10, seed = NULL, silent = FALSE, time = TRUE){
  st0 <- Sys.time()
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  #### Some calculation for kinship matrix ####
  if(is.null(K)){
    K <- A.mat(x)
  }
  
  if(is.null(n.lines)){
    n.lines <- nrow(K)
  }
  
  #### Some modification of phenotype data ####
  Y <- as.matrix(Y)
  num.pheno <- ncol(Y)
  
  if(CV){
    #### n-fold cross validation ####
    id <- sample(1:n.lines %% n.fold)
    id[id == 0] <- n.fold
    
    ids <- rep(id, env.num)
    if(pred.method == "OCR2"){
      Y.train.no <- (1:nrow(Y))[-Y.test.no]
      
      Y2 <- Y
      Y[Y.train.no, ncol(Y)] <- Y2[Y.test.no, ncol(Y)]
      Y[Y.test.no, ncol(Y)] <- Y2[Y.train.no, ncol(Y)]
      
      y.obs <- (Y[, ncol(Y)])[Y.test.no]
      y.pred <- rep(NA, length(Y.test.no))
    }else{
      if(pred.method == "OCR"){
        Y.train.no <- (1:nrow(Y))[-Y.test.no]
        y.obs <- (Y[, ncol(Y)])[Y.test.no]
        y.pred <- rep(NA, length(Y.test.no))
      }else{
        y.obs <- Y[, ncol(Y)]
        y.pred <- rep(NA, nrow(Y))
      }
    }
    m.ress <- NULL
    
    #### Prepare the training and test data-set ####
    if(num.pheno >= 2){
      if(is.null(XF)){
        covariate <- cbind(rep(1, nrow(Y)), Y[, -ncol(Y)])
      }else{
        covariate <- cbind(rep(1, nrow(Y)), Y[, -ncol(Y)], XF)
      }
    }else{
      if(is.null(XF)){
        covariate <- cbind(rep(1, nrow(Y)))
      }else{
        covariate <- cbind(rep(1, nrow(Y)), XF)
      }
    }
    
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
      
      
      y.test <- Y[, ncol(Y)]
      if((pred.method == "OCR") | (pred.method == "OCR2")){
        y.NA.no <- sort(c(Y.train.no[ids == i], Y.test.no))
        y.test[y.NA.no] <- NA
      }else{
        y.test[ids == i] <- NA
      }
      
      
      #### The start of prediction by mixed.solve ####
      m.res <- mixed.solve(y = y.test, Z = Z, K = K, X = covariate, method = method,
                           bounds = bounds, SE = SE, return.Hinv = return.Hinv)
      
      if(num.pheno >= 2 | (!is.null(XF))){
        if((pred.method == "OCR") | (pred.method == "OCR2")){
          y.pred[ids == i] <- (covariate[Y.test.no[ids == i], ] %*% m.res$beta)[, 1] + 
            rep(m.res$u[id == i], env.num)
        }else{
          y.pred[ids == i] <- (covariate[ids == i, ] %*% m.res$beta)[, 1] + rep(m.res$u[id == i], env.num)
        }
      }else{
        if((pred.method == "OCR") | (pred.method == "OCR2")){
          y.pred[ids == i] <- rep(m.res$beta, sum(ids == i)) + rep(m.res$u[id == i], env.num)
        }else{
          y.pred[ids == i] <- rep(m.res$beta, sum(ids == i)) + rep(m.res$u[id == i], env.num)
        }
        
      }
      m.ress <- c(m.ress, list(m.res))
      end <- Sys.time()
    }
  }else{
    #### When you don't perform cross-validation ####
    if(is.null(Y.test.no)){
      stop("You should assign the test data-set no. to Y.test.no when you don't perform cross-validation!")
    }
    
    y.obs <- Y[Y.test.no, ncol(Y)]
    y.pred <- rep(NA, length(Y.test.no))
    
    y.test <- y.obs
    y.test[Y.test.no] <- NA
    
    if(is.null(XF)){ 
      covariate <- cbind(rep(1, nrow(Y)), Y[, -ncol(Y)])
    }else{
      covariate <- cbind(rep(1, nrow(Y)), Y[, -ncol(Y)], XF)
    }
    
    
    #### The start of prediction by MTM ####
    m.res <- mixed.solve(y = y.test, Z = Z, K = K, X = covariate, method = method,
                         bounds = bounds, SE = SE, return.Hinv = return.Hinv)
    
    if(num.pheno >= 2 | (!is.null(XF))){
      y.pred <- (covariate[Y.test.no, ] %*% m.res$beta)[,1] + m.res$u[Y.test.no]
    }else{
      y.pred <- rep(m.res$beta, sum(id == i)) + m.res$u[Y.test.no]
    }
    m.ress <- m.res
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
      if(pred.method == "OCR" || pred.method == "OCR2"){
        y.obs.now <- y.obs[((i - 1) * n.lines + 1):(i * n.lines)]
        y.pred.now <- y.pred[((i - 1) * n.lines + 1):(i * n.lines)]
      }else{
        y.obs.now <- y.obs[((i - 1) * nrow(Y) / 3 + 1):(i * nrow(Y) / 3)]
        y.pred.now <- y.pred[((i - 1) * nrow(Y) / 3 + 1):(i * nrow(Y) / 3)]
      }
      R2.now <- try(cor(y.obs.now, y.pred.now, use = "complete.obs") ^ 2, silent = TRUE)
      
      if(class(R2.now) == "try-error"){
        R2.now <- NA
      }
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
    
    R2 <- mean(R2s, na.rm = T)
    MSE <- mean(MSEs, na.rm = T)
    NDCG <- mean(NDCGs, na.rm = T)
  }
  
  
  
  end0 <- Sys.time()
  if(time){
    print(paste0("It took ", end0- st0, " ", attr(end0- st0, "units"), " for this prediction."))
  }
  
  return(list(m.res = m.ress, y.obs = y.obs, y.pred = y.pred, R2 = R2, MSE = MSE, NDCG = NDCG))
}




rrBLUP_prediction_iter <- function(Y, Y.test.no = NULL, XF = NULL, x = NULL, Z = NULL, K = NULL, n.iter = 10, fitting.scale = "R2",
                                   method = "REML", bounds = c(1e-09, 1e+09), SE = TRUE, return.Hinv = TRUE, NDCG.k = 10, 
                                   n.lines = NULL, env.num = 1, pred.method = "OCR2",
                                   CV = TRUE, n.fold = 10, seeds = NULL, silent = TRUE, time = TRUE, saveAt = NULL){
  if(!is.null(seeds)){
    if(length(seeds) != n.iter){
      stop("Please set so that the length of seeds argument is same as n.iter argument!")
    }
  }
  R2s <- MSEs <- NDCGs <- rep(NA, n.iter)
  for(i in 1:n.iter){
    res.now <- rrBLUP_prediction(Y = Y, Y.test.no = Y.test.no, XF = XF, x = x, Z = Z, K = K,
                                 method = method, bounds = bounds, SE = SE, return.Hinv = return.Hinv, NDCG.k = NDCG.k,
                                 n.lines = n.lines, env.num = env.num, pred.method = pred.method, 
                                 CV = CV, n.fold = n.fold, seed = seeds[i], silent = silent, time = time)
    
    
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
