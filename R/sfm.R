sfm <- function(formula, 
                model_name    = c("NHN","NHN_Z","NE","NE_Z","NR","THT","NTN","NG","NNAK"),
                data, 
                maxit.bobyqa  = 10000,
                maxit.psoptim = 1000,
                maxit.optim   = 1000, 
                REPORT        = 1, 
                trace         = 2,
                pgtol         = 0, 
                start_val     = FALSE,
                PSopt         = FALSE,
                optHessian    = TRUE,
                inefdec       = TRUE,
                upper         = NA,
                Method        = "L-BFGS-B",
                eta           = 0.01,   
                alpha         = 0.2,
                verbose       = FALSE,
                rand.psoptim  = NULL){

call <- match.call()
model_name    <- match.arg(model_name)   

DR1 <- data_proc(formula, data, model_name, individual = NULL, inefdec)

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

Start_Cs <- start_cs( formula_x ,data_orig, x_vars_vec, intercept, model_name, n_x_vars, start_val,n_z_vars,z_vars) 

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

if(model_name %in% c("NHN","NE","NR","NG","NNAK","THT","NTN","NHN_Z","NE_Z") ){
like.fn = function(x){
      
if(model_name %in% c("NHN","NE","NR","NG","NNAK")){x_x_vec <- x[3:as.numeric(n_x_vars + 2)]}
if(model_name %in% c("THT","NTN")){                                                x_x_vec <- x[4:as.numeric(n_x_vars + 3)]}

if(model_name %in% c("NE_Z","NHN_Z")){data_z_vars <- as.matrix(data.frame(subset(data,select = z_vars)))
                                      x_x_vec     <- x[2:as.numeric(n_x_vars+1)]
                                      z_z_vec     <- x[as.numeric(n_x_vars+2):as.numeric(length(start_v))]}

      eps     <- (inefdec_n*(Y  - as.matrix(data_i_vars)%*%x_x_vec))
      
      if(model_name=="NHN_Z"){
        sigma_u_fun    <- exp(as.matrix(data_z_vars)%*%z_z_vec)
        sigma_v_fun    <- x[1]
        sigma_fun      <- sqrt(sigma_v_fun^2 + sigma_u_fun^2)
        lamb_fun       <- sigma_u_fun/sigma_v_fun
        like           <-       log(pmax(   (2/sigma_fun)  * 
                                              dnorm(eps/sigma_fun)*  
                                              pnorm(-eps*lamb_fun/sigma_fun), .Machine$double.xmin) )}
      
      if(model_name=="NE_Z"){
        sigma_u_fun    <- exp(data_z_vars%*%z_z_vec)
        sigv           <- x[1]
        l1             <- log(1/sigma_u_fun)
        l2             <- pnorm( -(eps/sigv) - (sigv    /sigma_u_fun), log.p = TRUE)
        l3             <- (eps/sigma_u_fun) + (sigv^2 /  (2*sigma_u_fun^2) )
        like           <-  l1+l2+l3}
      
      if(model_name == "NHN"){
      like  <-      as.numeric(log(     pmax(   (2/x[2])    * 
                                              dnorm(eps/x[2]) *
                                              pnorm(-eps*x[1]/x[2])  ,  .Machine$double.xmin )    ))}
      
      if(model_name == "NE"){
      l1   <- log(1/x[2])
      l2   <- pnorm( -(eps/x[1]) - (x[1]    /x[2]), log.p = TRUE)
      l3   <- (eps/x[2]) + (x[1]^2 /  (2*x[2]^2)  )
      like <-  l1+l2+l3}
      
      if(model_name=="NR"){
        sigv           <- x[1]
        sigu           <- x[2]
        sigma          <- sqrt(2*sigv^2+sigu^2)
        z              <- (eps*sigu/sigv)/sigma
        like           <- (log(pmax(sigv,.Machine$double.xmin))- 2*log(pmax(sigma,.Machine$double.xmin))
                          - 1/2*(eps/sigv)^2 + 1/2*z^2 + log(pmax(sqrt(2/pi)*exp(-1/2*z**2)
                          - z*(1-erf(z/sqrt(2))),.Machine$double.xmin)))}
      
      if(model_name=="THT"){
        sig_u   <- x[1]
        sig_v   <- x[2]
        a       <- x[3]
        lamb    <- sig_u/sig_v
        sig     <- sqrt(sig_v^2 + sig_u^2)
        like    <- as.numeric(log(pmax(2*dt(eps, df=a)*
                   pt((-eps*lamb/sig)*sqrt((a+1)/(sig^{-2}*eps^2 + a))  ,df=a+1)  , .Machine$double.xmin))) }
      
      if(model_name=="NTN"){
        lam  <- x[1]
        sig  <- x[2]
        mu   <- x[3]
        
        l1   <- -log(sig^2)/2
        l2   <- -log(2*pi)/2
        l3   <- -(1/(2*sig^2))*(-eps-mu)^2  
        l4   <-  pnorm(((mu/lam)-eps*lam)/sig,   log.p=TRUE)  
        l5   <- -pnorm((mu/sig)*sqrt(1+lam^(-2)),log.p=TRUE)  
        like <- l1 + l2 + l3 + l4 + l5}

      if(model_name %in% c("NG","NNAK")){
        lnDv <- function(nu,z){
          ((nu/2)*log(2) + 0.5*log(pi) -z^2/4 + log(1/gamma((1-nu)/2)*hyperg_1F1(-nu/2,1/2,z^2/2)
          -sqrt(2)*z/gamma(-nu/2)*hyperg_1F1((1-nu)/2,3/2,z^2/2)))}
        sig_v  <- x[1]
        sig_u  <- x[2]
        mu     <- x[3]
        if(model_name=="NG"){
          like   <- ((mu-1)*log(sig_v) - 1/2*log(2) - 1/2*log(pi) - mu*log(sig_u)
                    - 1/2*(eps/sig_v)^2 + 1/4*(eps/sig_v+sig_v/sig_u)^2
                    + lnDv(-mu,eps/sig_v+sig_v/sig_u))}
        if(model_name=="NNAK"){
          sigma <- sqrt(2*mu*sig_v^2+sig_u^2)
          like   <- (lgamma(2*mu) - lgamma(mu) + 1/2*log(2) - 1/2*log(pi) + mu*log(mu)
                    + (2*mu-1)*log(sig_v) - 2*mu*log(sigma) - 1/2*(eps/sig_v)^2
                    + 1/4*((eps*sig_u/sig_v)/sigma)^2 + lnDv(-2*mu,(eps*sig_u/sig_v)/sigma))}}
      
      like[like==-Inf]         <-  -sqrt(.Machine$double.xmax/length(like))
      like[like== Inf]         <-  -sqrt(.Machine$double.xmax/length(like))
      like[is.nan(like)]       <-  -sqrt(.Machine$double.xmax/length(like))
      
      return(-sum(like[is.finite(like)]) ) }  
  
Start.Time <- start.time()
   
Opt.Bobyqa  <- opt.bobyqa(fn=like.fn, start_v=start_v, lower.bobyqa=lower_bob, maxit.bobyqa=maxit.bobyqa, bob.TF=TRUE,verbose = verbose) 
start_v     <- Opt.Bobyqa$start_v
start_feval <- Opt.Bobyqa$start_feval
bob1        <- Opt.Bobyqa$bob1 

Lower.Start <- lower.start(Opt.Bobyqa$start_v, model_name, differ=1)

Opt.Psoptim <- opt.psoptim(fn=like.fn, start_v, lower.psoptim=Lower.Start$lower1, rand.psoptim = rand.psoptim,
                   upper.psoptim=Lower.Start$upper1, maxit.psoptim, psopt.TF=PSopt, rand.order = FALSE,verbose = verbose)  
start_v     <- Opt.Psoptim$start_v
start_feval <- Opt.Psoptim$start_feval
opt00       <- Opt.Psoptim$opt00

Lower.Start <- lower.start(start_v, model_name, differ=0.5)
Opt.Optim   <-  opt.optim(fn = like.fn, start_v = start_v, lower.optim =Lower.Start$lower1,
      upper.optim=Lower.Start$upper1, maxit.optim=maxit.optim, opt.TF=optHessian, method=Method, optHessian= TRUE,verbose = verbose)
start_v     <- Opt.Optim$start_v
start_feval <- Opt.Optim$start_feval
opt         <- Opt.Optim$opt

End.Time <- end.time(Start.Time)    

if(optHessian==FALSE & PSopt == FALSE){opt <- bob1
st_err  <- rep(NA,length(opt$par))}

if(optHessian==FALSE & PSopt == TRUE){opt <- opt00
st_err  <- rep(NA,length(opt$par))}

if(optHessian==TRUE){ st_err <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0)){rep(NA,length(opt$par))}else{suppressWarnings(sqrt(diag(solve(opt$hessian))))}}
try(t_val      <- opt$par/st_err , silent = T)
out[1,]        <- opt$par
try(out[2,]    <- st_err         , silent = T)
try(out[3,]    <- t_val          , silent = T)

## JLMS TE Measurements     
if(model_name %in% c("NHN")){
beta      <- opt$par[-c(1:2)]
lamb      <- opt$par[1]
sig       <- opt$par[2]
sig_u     <- (lamb*sig) / sqrt(1+lamb^2)
sig_v     <- sig_u/lamb
eps_hat   <- inefdec_n*(Y - rowSums(t(t(data_i_vars)*beta))) 
sig_star  <- sig_u*sig_v/sig
inner     <- (lamb*eps_hat)/sig
exp_u_hat <- ((1-pnorm((sig_u*sig_v/sig) + inner ))/ pmax((1-pnorm(inner)), .Machine$double.xmin))*
             exp((sig_u^2/sig^2)*(eps_hat + 0.5*sig_v^2))
exp_u_hat <- pmax(exp_u_hat, 0)
exp_u_hat <- pmin(exp_u_hat, 1)
## Median TE 
mu        <- 0
sigu2     <- sig_u^2
sigv2     <- sig_v^2
sig2      <- sigv2+sigu2
mu.star   <- (sigv2*mu-sigu2*(eps_hat))/sig2
sig.star  <- sig_u*sig_v/sqrt(sig2)

med_u_hat <- exp(-mu.star+sig.star*qnorm(0.5*pnorm(mu.star/sig.star)))
med_u_hat <- pmax(med_u_hat, 0)
med_u_hat <- pmin(med_u_hat, 1)}

if(model_name=="NR"){
beta      <- opt$par[-c(1:2)]
sig_v     <- opt$par[1]
sig_u     <- opt$par[2]
eps_hat   <- inefdec_n*(Y - rowSums(t(t(data_i_vars)*beta))) 
sigma     <- sqrt(2*sig_v^2 + sig_u^2)
z         <- (eps_hat*sig_u/sig_v)/sigma
exp_u_hat <- (exp(1/2*(z+sig_v*sig_u/sigma)^2-1/2*z^2)*
              (exp(-1/2*(z+sig_v*sig_u/sigma)^2)-sqrt(pi/2)*(z+sig_v*sig_u/sigma)*(1-erf(1/sqrt(2)*(z+sig_v*sig_u/sigma))))/
              (exp(-z^2/2)-sqrt(pi/2)*z*(1-erf(z/sqrt(2)))))
exp_u_hat <- pmax(exp_u_hat, 0)}

if(model_name=="NHN_Z"){
NX         <- n_x_vars + 1
NZ1        <- n_x_vars + 2
NZ2        <- n_x_vars + n_z_vars + 1 
beta       <- opt$par[c(2:NX)]
delta      <- opt$par[c(NZ1:NZ2)]
sig_v      <- opt$par[1]
sig_u      <- exp((as.matrix(as.matrix(data.frame(subset(data,select = z_vars)))))%*%delta)
lamb       <- sig_u/sig_v 
sig        <- sqrt(sig_u^2 + sig_v^2)
eps_hat    <- inefdec_n*(Y - rowSums(t(t(data_i_vars)*beta)))
sig_star   <- (sig_u*sig_v)  /sig
inner      <- (lamb*eps_hat) /sig
exp_u_hat  <- ((1-pnorm((sig_u*sig_v/sig) + inner ))/  pmax( (1-pnorm(inner)), .Machine$double.xmin)  )*exp( (sig_u^2/sig^2)*(eps_hat + 0.5*sig_v^2) )
exp_u_hat  <- pmax(exp_u_hat, 0)
exp_u_hat  <- pmin(exp_u_hat, 1)}

if(model_name %in% c("NG","NNAK")){
lnDv <- function(nu,z){
  ((nu/2)*log(2) + 0.5*log(pi) -z^2/4 + log(1/gamma((1-nu)/2)*hyperg_1F1(-nu/2,1/2,z^2/2)
  -sqrt(2)*z/gamma(-nu/2)*hyperg_1F1((1-nu)/2,3/2,z^2/2),))}
beta  <- opt$par[-c(1:3)]
sig_v <- opt$par[1]
sig_u <- opt$par[2]
mu    <- opt$par[3]
eps_hat    <- inefdec_n*(Y - rowSums(t(t(data_i_vars)*beta)))
if(model_name=="NG"){
  z         <- eps_hat/sig_v + sig_v/sig_u
  exp_u_hat <- exp(((z+sig_v)/2)^2)/exp((z/2)^2)*exp(lnDv(-mu,z+sig_v))/exp(lnDv(-mu,z))}
if(model_name=="NNAK"){
  sigma     <- sqrt(2*mu*sig_v^2 + sig_u^2)
  z         <- (eps_hat*sig_u/sig_v)/sigma
  exp_u_hat <- exp((z/2 + sig_v*sig_u/(2*sigma))^2)/exp((z/2)^2)*exp(lnDv(-2*mu,z+sig_v*sig_u/sigma))/exp(lnDv(-2*mu,z))}
exp_u_hat  <- pmax(exp_u_hat, 0)
exp_u_hat  <- pmin(exp_u_hat, 1)}

if(model_name %in%  c("NHN") ){
results <- list(t(out),c(opt),End.Time,start_v,model_name,formula, exp_u_hat, med_u_hat, 
                out["par",], out["st_err",], out["t-val",],call  )
class(results) <- "sfareg"
names(results)  <- c("out","opt","total_time","start_v","model_name","formula","exp_u_hat","med_u_hat",
                     "coefficients", "std.errors", "t.values","call")}

if(model_name %in%  c("NHN_Z","NR","NG","NNAK") ){
results <- list(t(out),c(opt),End.Time,start_v,model_name,formula, exp_u_hat,
                out["par",], out["st_err",], out["t-val",],call )
class(results) <- "sfareg"
names(results)  <- c("out","opt","total_time","start_v","model_name","formula","exp_u_hat",
                     "coefficients", "std.errors", "t.values","call")}

if(model_name %nin%  c("NHN","NHN_Z","NR","NG","NNAK") ){
results <- list(t(out),c(opt),End.Time,start_v,model_name,formula,
                  out["par",], out["st_err",], out["t-val",],call )
class(results) <- "sfareg"
names(results)  <- c("out","opt","total_time","start_v","model_name","formula",
                     "coefficients", "std.errors", "t.values","call")}
return(results)}

else {return(c("This is not a valid command"))}}



