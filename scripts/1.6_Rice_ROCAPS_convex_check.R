convex_check <- function(x, y, selected, return.all = FALSE){
  x.selected <- x[selected]
  y.selected <- y[selected]
  
  sort.x <- order(x.selected)
  x.sort <- x.selected[sort.x]
  y.sort <- y.selected[sort.x]
  sel.sort <- selected[sort.x]
  
  x.uniq <- unique(x.sort)
  y.uniq <- convex.uniq <- rep(NA, length(x.uniq))
  for(i in 1:length(x.uniq)){
    y.uniq[i] <- max(y.sort[x.sort == x.uniq[i]])
    convex.uniq[i] <- sel.sort[which(y.sort == max(y.sort[x.sort == x.uniq[i]]))]
  }
  
  slope.deltas <- 1

  while(sum(slope.deltas > 0) >= 1){
    slopes <- rep(NA, length(x.uniq) - 1)
    for(i in 1:(length(x.uniq) - 1)){
      slopes[i] <- (y.uniq[i + 1] - y.uniq[i]) / (x.uniq[i + 1] - x.uniq[i])
    }
    
    slope.deltas <- slope.ratios <- rep(NA, length(x.uniq) - 2)
    for(i in 1:(length(x.uniq) - 2)){
      slope.deltas[i] <- slopes[i + 1] - slopes[i]
      slope.ratios[i] <- slopes[i + 1]/slopes[i]
    }
    
    if(sum(slope.deltas > 0) >= 1){
      non.convex <- which(slope.deltas > 0) + 1
      convex.no <- convex.uniq[-non.convex]
    }else{
      convex.no <- convex.uniq
    }
    
    y.uniq <- y[convex.no]
    x.uniq <- x[convex.no]
    convex.uniq <- convex.no
  }
  
  
  if(!return.all){
    return(convex.no)
  }else{
    return(list(convex.no = convex.no, slopes = slopes, sl.deltas = slope.deltas, sl.ratio = slope.ratios))
  }
}
