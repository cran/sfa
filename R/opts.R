opt.bobyqa    <- function(fn, start_v, lower.bobyqa, maxit.bobyqa, bob.TF, rhobeg = NA, rhoend  =NA, verbose=verbose){
start_feval   <-  fn(start_v)
bob1          <- NULL 
if(isTRUE(bob.TF==TRUE)){  

bob1   <- bobyqa(par = start_v, 
                 fn = fn,
                 lower   = lower.bobyqa, 
                 control = list(iprint = if (verbose) 2 else 0, 
                                maxfun = maxit.bobyqa,
                                rhobeg = rhobeg,
                                rhoend = rhoend))

if(isTRUE(start_feval > bob1$fval )) {start_v <- bob1$par
start_feval   <-  fn(start_v)} }

results        <- list(start_v,start_feval,   bob1)   
names(results) <-  c("start_v","start_feval","bob1")  
return(results)}

opt.optim     <- function(fn, start_v, lower.optim, upper.optim, maxit.optim, opt.TF, method, optHessian, trace, verbose=verbose){
  start_feval   <-  fn(start_v)
  opt           <- NULL 
  if(isTRUE(opt.TF ==TRUE)){
  
  
    
    opt <- optim(par     = start_v, 
                 fn      = fn,
                 lower   = lower.optim, 
                 upper   = upper.optim,
                 hessian = optHessian, 
                 method  = method,
                 control = list(maxit   = maxit.optim, 
                                REPORT  = base::ceiling(maxit.optim/10), 
                                trace   = if(verbose) {1} else {0}))
    
    if(isTRUE(start_feval > opt$value )) {start_v <- opt$par
    start_feval   <-  fn(start_v)} }
  
  results        <- list(start_v,start_feval,   opt)   
  names(results) <-  c("start_v","start_feval","opt")  
  return(results)}

opt.psoptim   <- function(fn, start_v, lower.psoptim, upper.psoptim=NA, maxit.psoptim=maxit.psoptim, 
                          psopt.TF, rand.order = TRUE, verbose=verbose, rand.psoptim=rand.psoptim){
  start_feval <-  fn(start_v)
  opt00       <- NULL 
  report      <- base::ifelse(verbose, base::ceiling(maxit.psoptim/10), 0)
  if(isTRUE(psopt.TF ==TRUE)){  
    
    if(!is.null(rand.psoptim)){
       set.seed(rand.psoptim)}
    
    opt00 <- psoptim(par     = start_v, 
                     fn      = fn,
                     lower   = lower.psoptim , 
                     upper   = upper.psoptim,
                     control = list(trace          = if(verbose){1} else{0},
                                    REPORT         = report,
                                    trace.stats    = if(verbose){TRUE} else{FALSE},
                                    maxit          = maxit.psoptim,
                                    rand.order     = rand.order))
    
    if(isTRUE(start_feval > opt00$value )) {start_v <- opt00$par
    start_feval   <-  fn(start_v)}  }
  
  results        <- list(start_v,start_feval,   opt00)   
  names(results) <-  c("start_v","start_feval","opt00")  
  return(results)}
