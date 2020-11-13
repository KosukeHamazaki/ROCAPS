Acquisition_func <- function(x, y, selected, m, sd){
  x.selected <- x[selected]
  y.selected <- y[selected]
  
  sort.x <- order(x.selected)
  x.sort <- x.selected[sort.x]
  y.sort <- y.selected[sort.x]
  
  y.thres <- acquisition <- rep(NA, length(x))
  for(i in 1:(length(x.sort) - 1)){
    slope <- (y.sort[i + 1] - y.sort[i]) / (x.sort[i + 1] - x.sort[i])
    intercept <- y.sort[i] - slope * x.sort[i]
    
    y.thres[(x >= x.sort[i]) & (x < x.sort[i + 1])] <- x[(x >= x.sort[i]) & (x < x.sort[i + 1])] * slope + intercept
    acquisition[(x >= x.sort[i]) & (x < x.sort[i + 1])] <- pnorm(q = y.thres[(x >= x.sort[i]) & (x < x.sort[i + 1])],
                                                                        mean = m[(x >= x.sort[i]) & (x < x.sort[i + 1])],
                                                                        sd = sd[(x >= x.sort[i]) & (x < x.sort[i + 1])],
                                                                        lower.tail = FALSE)
    
  }
  y.thres[length(x)] <- y.sort[length(x.selected)]
  acquisition[length(x)] <-  pnorm(q = y.thres[length(x)], mean = m[length(x)], sd = sd[length(x)], lower.tail = FALSE)
  return(list(Acquire = acquisition, y.thres = y.thres))
}
