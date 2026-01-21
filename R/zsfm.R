zsfm <- function(formula, 
                model_name    = c("ZISF","ZISF_Z"),
                data, 
                maxit.bobyqa  = 10000,
                maxit.psoptim = 1000,
                maxit.optim   = 1000, 
                REPORT        = 1, 
                trace         = 0,
                pgtol         = 0, 
                start_val     = FALSE,
                PSopt         = FALSE,
                optHessian    = TRUE,
                inefdec       = TRUE,
                upper         = NA,
                Method        = "L-BFGS-B",
                logit         = TRUE,
                verbose       = FALSE,
                rand.psoptim  = NULL){
  
call <- match.call()
model_name    <- match.arg(model_name) 

DR1 <- data_proc(formula,   data, model_name, individual = NULL, inefdec)

data_orig    <- DR1$data_orig
form_parts   <- DR1$form_parts
formula_x    <- DR1$formula_x
y_var        <- DR1$y_var
model_name   <- DR1$model_name
data_x       <- DR1$data_x
intercept    <- DR1$intercept
inefdec_n    <- DR1$inefdec_n
inefdec_TF   <- DR1$inefdec_TF
x_vars_vec   <- DR1$x_vars_vec
n_x_vars     <- DR1$n_x_vars
x_vars       <- DR1$x_vars
x_x_vec      <- DR1$x_x_vec
fancy_vars   <- DR1$fancy_vars
fancy_vars_z <- DR1$fancy_vars_z
n_z_vars     <- DR1$n_z_vars
N            <- DR1$N
data_z       <- DR1$data_z
if(length(unlist(form_parts))>3){  
formula_z    <- DR1$formula_z
intercept_z  <- DR1$intercept_z
n_z_vars     <- DR1$n_z_vars
z_vars       <- DR1$z_vars
z_vars_vec   <- DR1$z_vars_vec
z_z_vec      <- DR1$z_z_vec}


Start_Cs     <- start_cs( formula_x ,data_orig, x_vars_vec, intercept, model_name, n_x_vars, start_val,n_z_vars,z_vars) 

beta_0         <- Start_Cs$beta_0
beta_0_st      <- Start_Cs$beta_0_st
beta_hat       <- Start_Cs$beta_hat
epsilon_hat    <- Start_Cs$epsilon_hat
lambda         <- Start_Cs$lambda
lower_bob      <- Start_Cs$lower_bob
mu             <- Start_Cs$mu
out            <- Start_Cs$out
plm_lm         <- Start_Cs$plm_lm
sigma          <- Start_Cs$sigma
sigma_u        <- Start_Cs$sigma_u
sigma_v        <- Start_Cs$sigma_v
start_v        <- Start_Cs$start_v
start_v_ne     <- Start_Cs$start_v_ne
start_v_ng     <- Start_Cs$start_v_ng
start_v_nhn    <- Start_Cs$start_v_nhn
start_v_nnak   <- Start_Cs$start_v_nnak
start_v_nr     <- Start_Cs$start_v_nr
start_v_ntn    <- Start_Cs$start_v_ntn
start_v_t      <- Start_Cs$start_v_t

DR2 <- data_proc2(data, data_x, fancy_vars, fancy_vars_z, data_z, y_var, x_vars_vec, halton_num=NA, individual=NA, N, model_name, rand.gtre=NULL)

data         <- DR2$data
Y            <- DR2$Y
data_i_vars  <- DR2$data_i_vars


if(model_name %in% c("ZISF","ZISF_Z") ){
like.fn = function(x){
      
      if(model_name %in% c("ZISF")){  x_x_vec <- x[4:as.numeric(n_x_vars+3)]}
      if(model_name %in% c("ZISF_Z")){x_x_vec <- x[3:as.numeric(n_x_vars+2)]}
        
      eps     <- (inefdec_n*(Y  - as.matrix(data_i_vars)%*%x_x_vec))
      
if(model_name == "ZISF"){
  	gamma  <- x[1]
  	prob   <- exp(-abs(gamma)) 
  	sigvsq <- x[2]^2
  	sigusq <- x[3]^2
  	sigv   <- sqrt(sigvsq)
  	sigu   <- sqrt(sigusq)

  	lambda <- sigu/sigv
  	sigsq  <- sigvsq+sigusq
  	sig    <- sqrt(sigsq)           

  	f1     <- -0.5*log(2*pi*sigvsq)-(0.5/sigvsq)*eps^2
  	f2     <- log(2/sig)+log(dnorm(eps/sig))+log(pnorm(eps*lambda/sig))
  	f      <- prob*exp(f1)+(1-prob)*exp(f2)

  	like   <- log(f+1e-10) }
      
if(model_name == "ZISF_Z"){
    gamma <- x[(n_x_vars+3):(n_x_vars+2+n_z_vars)]  ## lets put gammas last 
    
    if(logit){ prob <- exp(  data_z%*%gamma)/(1+exp(  data_z%*%gamma))}
    if(!logit){prob <- pnorm(data_z%*%gamma)/(1+pnorm(data_z%*%gamma))}
    
    sigvsq <- x[1]^2
    sigusq <- x[2]^2
    sigv   <- sqrt(sigvsq)
    sigu   <- sqrt(sigusq)
    
    lambda <- sigu/sigv
    sigsq  <- sigvsq+sigusq
    sig    <- sqrt(sigsq)           
      
    f1     <- -0.5*log(2*pi*sigvsq)-(0.5/sigvsq)*eps^2
    f2     <- log(2/sig)+log(dnorm(eps/sig))+log(pnorm(eps*lambda/sig))
    f      <- prob*exp(f1)+(1-prob)*exp(f2)
        
    like   <- log(f+1e-10) }      


like[like==-Inf]         <-  -sqrt(.Machine$double.xmax/length(like))
like[like== Inf]         <-  -sqrt(.Machine$double.xmax/length(like))
like[is.nan(like)]       <-  -sqrt(.Machine$double.xmax/length(like))
      
return(-sum(like[is.finite(like)]))}  
  
Start.Time <- start.time()
   
Opt.Bobyqa <- opt.bobyqa(fn=like.fn, start_v=start_v, lower.bobyqa=lower_bob, maxit.bobyqa=maxit.bobyqa, bob.TF=TRUE,verbose = verbose) 
start_v     <- Opt.Bobyqa$start_v
start_feval <- Opt.Bobyqa$start_feval
bob1        <- Opt.Bobyqa$bob1 

Lower.Start <- lower.start(start_v, model_name, differ=1)

Opt.Psoptim <- opt.psoptim(fn=like.fn, start_v, lower.psoptim=Lower.Start$lower1, rand.psoptim = rand.psoptim,
               upper.psoptim=Lower.Start$upper1, maxit.psoptim, psopt.TF=PSopt, rand.order = FALSE,verbose = verbose)  
start_v     <- Opt.Psoptim$start_v
start_feval <- Opt.Psoptim$start_feval
opt00       <- Opt.Psoptim$opt00

Lower.Start <- lower.start(start_v, model_name, differ=0.5)
Opt.Optim   <- opt.optim(fn = like.fn, start_v = start_v, lower.optim =Lower.Start$lower1,
          upper.optim=Lower.Start$upper1, maxit.optim=maxit.optim, opt.TF=optHessian, method=Method, optHessian= TRUE,verbose = verbose)
start_v     <- Opt.Optim$start_v
start_feval <- Opt.Optim$start_feval
opt         <- Opt.Optim$opt

End.Time <- end.time(Start.Time) 
if(optHessian==FALSE & PSopt == FALSE){opt <- bob1
st_err  <- rep(NA,length(opt$par))}

if(optHessian==FALSE & PSopt == TRUE){opt <- opt00
st_err  <- rep(NA,length(opt$par))}

if(optHessian==TRUE){ st_err  <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0 ) ){ rep(NA,length(opt$par)) }   else{suppressWarnings(sqrt(diag(solve(opt$hessian))))}}
t_val      <- opt$par/st_err
out[1,]    <- opt$par
out[2,]    <- st_err
out[3,]    <- t_val

## JLMS
if(model_name %in% c("ZISF","ZISF_Z")){

  if(is.na(n_z_vars)==TRUE){
  beta <- opt$par[-c(1:3)]
  z <- 1
  gamma <- opt$par[1]
  prob  <- exp(-gamma)
  
  sigvsq <- opt$par[2]^2
  sigusq <- opt$par[3]^2 }
  
  if(is.na(n_z_vars)==FALSE){
    beta  <- opt$par[3:sum(n_x_vars,2) ]
    gamma <- opt$par[(n_x_vars+3):(n_x_vars+2+n_z_vars)]  ## lets put gammas last 
    
    if(logit){prob  <- exp(data_z%*%gamma)/(1+exp(data_z%*%gamma))}
    if(!logit){prob <- pnorm(data_z%*%gamma)/(1+pnorm(data_z%*%gamma))}
    
    sigvsq <- opt$par[1]^2
    sigusq <- opt$par[2]^2 }
    
  eps  <- (inefdec_n*(Y  - as.matrix(data_i_vars)%*%beta))
  sigv <- sqrt(sigvsq)
  sigu <- sqrt(sigusq)
  
  ## Reparametrize the log-likelihood function
  lambda      <- sigu/sigv
  sigsq       <- sigvsq+sigusq
  sig         <- sqrt(sigsq)           
  
  ## Now the likelihood function
  f1          <- -0.5*log(2*pi*sigvsq)-(0.5/sigvsq)*eps^2
  f2          <- log(2/sig)+log(dnorm(eps/sig))+log(pnorm(eps*lambda/sig))
  f           <- prob*exp(f1)+(1-prob)*exp(f2)

  post.prob   <- prob*exp(f1)/f
  mustar      <- eps*sigusq/sigsq
  sigstarsq   <- sigusq*sigvsq/(sigusq+sigvsq)
  sigstar     <- sqrt(sigstarsq)
  zz          <- mustar/sigstar
  jlms        <- mustar+sigstar*dnorm(zz)/pnorm(zz) }

  if(model_name %in%  c("ZISF","ZISF_Z") ){
    results <- list(t(out),c(opt),End.Time,start_v,model_name,formula, jlms,post.prob,
                      out["par",], out["st_err",], out["t-val",],call)
    class(results)  <- "sfareg"
    names(results)  <- c("out","opt","total_time","start_v","model_name","formula","jlms","post.prob",
                           "coefficients", "std.errors", "t.values","call")}
    return(results)}
else {return(c("This is not a valid command"))}}



