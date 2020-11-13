require(BGLR)

GP <- function(y, K = NULL, X = NULL, error = FALSE, method = "RKHS", nIter = 6000, burnIn = 1000,
               thin = 5, saveAt = "", verbose = FALSE){
  if(method == "RKHS"){
    if(length(y) != nrow(K)){
      stop("The length of y should be equal to the dimension of K!")
    }
    
    if(error){
      E <- diag(length(y))
      ETA <- list(K = list(K = K, model = "RKHS"), E = list(K = E, model = "RKHS"))
    }else{
      ETA <- list(K = list(K = K, model = "RKHS"))
    }
  }else{
    if(method == "BRR"){
      ETA <- list(list(X = X, model = "BRR"))
    }
  }
  
  
  fm <- BGLR(y = y, ETA = ETA, nIter = nIter, burnIn = burnIn,
             saveAt = saveAt, thin = thin, verbose = verbose)
  
  
  return(list(mean = fm$yHat, sd = fm$SD.yHat))
}
