require(randomForest)


RF_prediction <- function(Y, Y.test.no = NULL, XF = NULL, x = NULL, Z = NULL, K = NULL,
                          n.lines = NULL, env.num = 1, pred.method = "OCR2",
                          pca.option = "prop", pca.prop = 0.5, pca.num = 10, NDCG.k = 10,
                          CV = TRUE, n.fold = 10, seed = NULL, silent = FALSE, time = TRUE){
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
  
  
  #### Prepare discriptional variable (using PCA) ####
  pca.res <- prcomp(K)
  
  if(pca.option == "prop"){
    pca.var.sum <- cumsum((pca.res$sdev / sum(pca.res$sdev)))
    pca.num <- min(which(pca.var.sum >= pca.prop))
  }
  
  pca.scores <- tcrossprod(Z %*% pca.res$x, Z)[, 1:pca.num] 
  
  if(num.pheno >= 2){
    if(is.null(XF)){
      setumei <- cbind(rep(1, nrow(Y)), Y[, -ncol(Y)], pca.scores)
    }else{
      setumei <- cbind(rep(1, nrow(Y)), Y[, -ncol(Y)], XF, pca.scores)
    }
  }else{
    if(is.null(XF)){
      setumei <- cbind(rep(1, nrow(Y)), pca.scores)
    }else{
      setumei <- cbind(rep(1, nrow(Y)), XF, pca.scores)
    }
  }
  
  if(CV){
    #### n-fold cros-svalidation ####
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
    
    RF.ress <- NULL
    
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
      
      #### Prepare the training and test data-set ####
      if((pred.method == "OCR") | (pred.method == "OCR2")){
        test.no <- sort(c(Y.train.no[ids == i], Y.test.no))
        x.train <- setumei[-test.no, ]
        x.test <- setumei[Y.test.no[ids == i], ]
        y.train <- Y[-test.no, ncol(Y)]
      }else{
        x.train <- setumei[ids != i, ]
        x.test <- setumei[ids == i, ]
        y.train <- Y[ids != i, ncol(Y)]
      }
      
      isna <- is.na(y.train)
      y.train <- y.train[!isna]
      x.train <- x.train[!isna, ]
      
      
      #### The start of prediction by Random Forest ####
      RF.res <- randomForest(x = x.train, y = y.train,  xtest = x.test,
                             ytest = NULL, ntree=500, importance = TRUE)
      
      y.pred[ids == i] <- RF.res$test$predicted
      
      RF.ress <- c(RF.ress, list(RF.res))
      end <- Sys.time()
    }
  }else{
    #### When you don't perform cross-validation ####
    if(is.null(Y.test.no)){
      stop("You should assign the test data-set no. to Y.test.no when you don't perform cross-validation!")
    }
    
    #### Prepare the training and test data-set ####
    x.train <- setumei[-Y.test.no, ]
    x.test <- setumei[Y.test.no, ]
    
    y.train <- Y[-Y.test.no, ncol(Y)]
    y.obs <- Y[Y.test.no, ncol(Y)]
    
    
    #### The start of prediction by Random Forest ####
    RF.res <- randomForest(x = x.train, y = y.train, xtest = x.test,
                           ytest = NULL, ntree=500, importance = TRUE)
    
    y.pred <- RF.res$test$predicted
    
    RF.ress <- RF.res
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
    print(paste0("It took ", end0- st0, " ", attr(end0- st0, "units"), " for this prediction."))
  }
  
  return(list(RF.res = RF.ress, y.obs = y.obs, y.pred = y.pred, R2 = R2, MSE = MSE, NDCG = NDCG))
}




RF_prediction_iter <- function(Y, Y.test.no = NULL, XF = NULL, x = NULL, Z = NULL, K = NULL, n.iter = 10, fitting.scale = "R2",
                               pca.option = "prop", pca.prop = 0.5, pca.num = 10, NDCG.k = 10,
                               n.lines = NULL, env.num = 1, pred.method = "OCR2",
                               CV = TRUE, n.fold = 10, seeds = NULL, silent = TRUE, time = TRUE, saveAt = NULL){
  if(!is.null(seeds)){
    if(length(seeds) != n.iter){
      stop("Please set so that the length of seeds argument is same as n.iter argument!")
    }
  }
  R2s <- MSEs <- NDCGs <- rep(NA, n.iter)
  for(i in 1:n.iter){
    res.now <- RF_prediction(Y = Y, Y.test.no = Y.test.no, XF = XF, x = x, Z = Z, K = K,
                             pca.option = pca.option, pca.prop = pca.prop, pca.num = pca.num,
                             n.lines = n.lines, env.num = env.num, pred.method = pred.method, NDCG.k = NDCG.k,
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
