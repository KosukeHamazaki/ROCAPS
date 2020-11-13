require(foreach)


stopImplicitCluster2 <- function() {
  options <- doParallel:::.options
  if (exists(".revoDoParCluster", where = options) &&
      !is.null(get(".revoDoParCluster", envir = options))) {
    stopCluster(get(".revoDoParCluster", envir = options))
    remove(".revoDoParCluster", envir = options)
  }
}

pforeach <- function(..., .combine=c, .parallel=TRUE, .debug=!.parallel, .cores, .seed=NULL, .export, .packages) {
  if(!require(doParallel)) stop("install.packages('doParallel')")
  if(!.parallel || .debug) {
    foreach::registerDoSEQ()
    if(!is.null(.seed)) set.seed(.seed)
  } else {
    if(missing(.cores)) .cores=parallel::detectCores()
    else if(.cores <= 0) .cores=parallel::detectCores() + .cores
    doParallel::registerDoParallel(.cores)
  }
  if(missing(.export)) .export=ls(parent.frame(1000))
  if(missing(.packages)) .packages=loadedNamespaces()
  return(function(expr) {
    expr <- substitute(expr)
    on.exit(stopImplicitCluster2())
    `%doop%` <- `%dopar%`
    if(!is.null(.seed)) {
      if(!require(doRNG)) stop("install.packages('doRNG')")
      set.seed(.seed)
      `%doop%` <- `%dorng%`
    }
    foreach(..., .combine=.combine, .export=.export, .packages=.packages) %doop% eval(expr)
  })
}