data_scaled_imputation <- function(Y, XF = NULL, Z, imp = TRUE, rm.NA = TRUE, scale.target = TRUE){
  if(nrow(Y) != nrow(Z)){
    stop("The number of lines doesn't match between Y and Z!")
  }

  Y.x <- Y[, -ncol(Y)]
  Y.y <- Y[, ncol(Y)] 
  
  Y.x.scaled <- scale(Y.x)
  if(scale.target){
    Y.y <- scale(Y.y)
  }

  
  
  if(imp){
    Y.x.scaled[is.na(Y.x.scaled)] <- 0
  }
  
  na.no <- which(is.na(Y.y))
    
  if(any(is.na(Y.y)) & rm.NA){
    Y.x.scaled2 <- Y.x.scaled[-na.no, ]
    Y.y2 <- Y.y[-na.no]
    Z2 <- Z[-na.no, ]
    
    if(!is.null(XF)){
      XF2 <- XF[-na.no, ]
    }else{
      XF2 <- XF
      }
  }else{
    Y.x.scaled2 <- Y.x.scaled
    Y.y2 <- Y.y
    Z2 <- Z
    XF2 <- XF
  }
  
  Y2 <- cbind(Y.x.scaled2, Y.y2)
  
  return(list(Y = Y2, XF = XF2, Z = Z2, na.no = na.no))
}


