start_cs    <- function(formula_x ,data_orig, x_vars_vec, intercept, model_name, n_x_vars, start_val,n_z_vars,z_vars){
plm_lm      <- lm(formula_x ,data_orig)  
beta_hat    <- if(isTRUE(intercept==0)) {plm_lm$coefficients[x_vars_vec]} else{plm_lm$coefficients[x_vars_vec][-1]}
epsilon_hat <- plm_lm$residuals
beta_0_st   <- if(isTRUE(intercept==0)) {NA} else{plm_lm$coefficients[c(1)]}
sigma_u     <- 0.1 
sigma_v     <- 0.1 
mu          <- 0.1 
beta_0      <- beta_0_st  
lambda      <- sigma_u/sigma_v
sigma       <- sqrt(sigma_u^2 + sigma_v^2)
gam         <- 1 

start_v_ntn    <- if(is.na(beta_0_st)) {unname(c(lambda,sigma,mu,beta_hat))}      else{unname(c(lambda,sigma,mu,beta_0,beta_hat)) }
start_v_t      <- if(is.na(beta_0_st)) {unname(c(sigma_u,sigma_v,1,beta_hat))}    else{unname(c(sigma_u,sigma_v,1,beta_0,beta_hat)) }
start_v_ng     <- if(is.na(beta_0_st)) {unname(c(sigma_v,sigma_u,1,beta_hat))}    else{unname(c(sigma_v,sigma_u,1,beta_0,beta_hat)) }
start_v_nnak   <- if(is.na(beta_0_st)) {unname(c(sigma_v,sigma_u,0.5,beta_hat))}  else{unname(c(sigma_v,sigma_u,0.5,beta_0,beta_hat)) }
start_v_ne     <- if(is.na(beta_0_st)) {unname(c(sigma_v,sigma_u,beta_hat))}      else{unname(c(sigma_v,sigma_u,beta_0,beta_hat)) }
start_v_nhn    <- if(is.na(beta_0_st)) {unname(c(lambda,sigma,beta_hat))}         else{unname(c(lambda,sigma,beta_0,beta_hat)) }
start_v_nr     <- if(is.na(beta_0_st)) {unname(c(sigma_v,sigma_u,beta_hat))}      else{unname(c(sigma_v,sigma_u,beta_0,beta_hat)) }
start_v_zisf   <- if(is.na(beta_0_st)) {unname(c(gam,sigma_v,sigma_u,beta_hat))}  else{unname(c(gam,sigma_v,sigma_u,beta_0,beta_hat)) }

if(model_name %in% c("NHN_Z","NE_Z")){
  start_v_nez  <- if(is.na(beta_0_st)) {unname(c(sigma_v,beta_hat,rep(mu,n_z_vars)))}else{unname(c(sigma_v,beta_0,beta_hat,rep(mu,n_z_vars))) }
  start_v        <- start_v_nez
  out            <- matrix(0,nrow = 3,ncol = length(start_v))
  colnames(out)  <- c("sigma_v",c(names(plm_lm$coefficients)),z_vars)
  lower_bob      <- c(.Machine$double.eps , rep(-.Machine$double.xmax^.1,length(start_v[-c(1)])))}
if(model_name %in% c("ZISF")){
  start_v        <- start_v_zisf
  out            <- matrix(0,nrow = 3,ncol = length(start_v))
  colnames(out)  <- c("gamma","sigv","sigu",c(names(plm_lm$coefficients)))
  lower_bob      <- c(-Inf, rep(.Machine$double.eps,2), rep(-Inf,n_x_vars) )}
if(model_name %in% c("ZISF_Z")){
start_v_zisfz  <- if(is.na(beta_0_st)) {unname(c(sigma_v,sigma_u,beta_hat,rep(mu,n_z_vars)))}else{unname(c(sigma_v,sigma_u,beta_0,beta_hat,rep(mu,n_z_vars))) }
  start_v        <- start_v_zisfz
  out            <- matrix(0,nrow = 3,ncol = length(start_v))
  colnames(out)  <- c("sigv","sigu",c(names(plm_lm$coefficients)),z_vars )
  lower_bob      <- c(rep(.Machine$double.eps,2), rep(-Inf,n_x_vars+n_z_vars) )}
if(model_name %in% c("NHN") ){
  start_v        <- start_v_nhn
  out            <- matrix(0,nrow = 3,ncol = length(start_v))
  colnames(out)  <- c("lambda","sigma",c(names(plm_lm$coefficients)))
  lower_bob      <- c(rep(.Machine$double.eps,2),rep(-Inf,n_x_vars))}
if(model_name=="NE"){
  start_v        <- start_v_ne
  out            <- matrix(0,nrow = 3,ncol = length(start_v))
  colnames(out)  <- c("sigv","sigu",c(names(plm_lm$coefficients)))
  lower_bob      <- c(rep(.Machine$double.eps,2),rep(-Inf,n_x_vars))}
if(model_name=="NR"){
  start_v        <- start_v_nr
  out            <- matrix(0,nrow = 3,ncol = length(start_v))
  colnames(out)  <- c("sigv","sigu",c(names(plm_lm$coefficients)))
  lower_bob      <- c(rep(.Machine$double.eps,2),rep(-Inf,n_x_vars))}
if(model_name=="THT"){
  start_v        <- start_v_t
  out            <- matrix(0,nrow = 3,ncol = length(start_v))
  colnames(out)  <- c("sigv","sigu","a",c(names(plm_lm$coefficients))) 
  lower_bob      <- c(rep(.Machine$double.eps,3),rep(-Inf,n_x_vars))}
if(model_name=="NG"){
  start_v        <- start_v_ng
  out            <- matrix(0,nrow = 3,ncol = length(start_v))
  colnames(out)  <- c("sigv","sigu","mu",c(names(plm_lm$coefficients))) 
  lower_bob      <- c(rep(.Machine$double.eps,3),rep(-Inf,n_x_vars))}
if(model_name=="NNAK"){
  start_v        <- start_v_nnak
  out            <- matrix(0,nrow = 3,ncol = length(start_v))
  colnames(out)  <- c("sigv","sigu","mu",c(names(plm_lm$coefficients))) 
  lower_bob      <- c(rep(.Machine$double.eps,3),rep(-Inf,n_x_vars))}
if(model_name=="NTN"){
  start_v        <- start_v_ntn
  out            <- matrix(0,nrow = 3,ncol = length(start_v))
  colnames(out)  <- c("lambda","sigma","mu",c(names(plm_lm$coefficients))) 
  lower_bob      <- c(rep(.Machine$double.eps,2),rep(-Inf,n_x_vars+1))}

if (isTRUE(is.numeric(start_val))) {start_v <- start_val} 

rownames(out)  <- c("par","st_err","t-val") 

results <- list(plm_lm, beta_hat, epsilon_hat, beta_0_st, sigma_u, sigma_v, mu,
             beta_0, lambda, sigma, start_v_ntn,start_v_ng, start_v_nnak,
             start_v_t, start_v_ne, start_v_nr, start_v_nhn, start_v, out,
             lower_bob)

names(results) <- c("plm_lm", "beta_hat", "epsilon_hat", "beta_0_st", "sigma_u", "sigma_v", "mu",
                 "beta_0", "lambda", "sigma", "start_v_ntn","start_v_ng", "start_v_nnak",
                 "start_v_t", "start_v_ne", "start_v_nr", "start_v_nhn", "start_v", "out",
                 "lower_bob")

return(results)
}

start_panel <- function(formula_x, data, model_name, start_val, intercept, x_vars_vec){
sfa_eps <- sfa_alp <- exp_eta <- exp_u <- sigma_v <- sigma_u <- 
sigma_r <- sigma_h <- beta_0 <- lambda <- sigma <- start_v <- out <- 
plm_gtre <- beta_hat <- alpha_hat <- epsilon_hat <- beta_0_st <- beta_se <- NULL

plm_tfe <- plm_fd <- NULL  

if (isFALSE(is.numeric(start_val))){
if(model_name %in% c("TFE","FD")){
plm_tfe     <- plm(formula_x,  data,effect = "individual",model = "within" )
plm_fd       <- plm(formula_x,  data,effect = "individual",model = "pooling")} else{
plm_gtre     <- plm(formula_x,  data,effect = "individual",model = "random" )
beta_hat     <- if(isTRUE(intercept==0)){plm(formula_x ,data,effect = "individual")$coefficients[c(x_vars_vec)]} else{plm_gtre$coefficients[-c(1)]}
alpha_hat    <- ranef(plm_gtre)
epsilon_hat  <- plm_gtre$residuals
beta_0_st    <- if(isTRUE(intercept==0)){NA}else{plm_gtre$coefficients[c(1)]}
beta_se      <- as.data.frame(summary(plm_gtre)[1])$coefficients.Std..Error
}}

if(isTRUE(is.numeric(start_val))){start_v <- start_val}
  
if(isFALSE(is.numeric(start_val)) &  model_name %in% c("TRE_Z","GTRE_Z","TRE","GTRE","GTRE_SEQ1","GTRE_SEQ2")){
  sfa_eps        <- pcs_c(Y  = as.numeric(epsilon_hat))[[1]]$par
  sfa_alp        <- pcs_c(Y  = as.numeric(alpha_hat  ))[[1]]$par
  exp_eta        <- sfa_alp[3]
  exp_u          <- sfa_eps[3]
  sigma_v        <- sqrt(sfa_eps[2]^2/ (1 + sfa_eps[1]^2))    
  sigma_u        <- sigma_v*sfa_eps[1]                        
  sigma_r        <- sqrt(sfa_alp[2]^2/ (1 + sfa_alp[1]^2))        
  sigma_h        <- sigma_r*sfa_alp[1]                        
  beta_0         <- beta_0_st + exp_u +exp_eta
  lambda         <- sigma_u/sigma_v
  sigma          <- sqrt(sigma_u^2 + sigma_v^2)
if(model_name %in% c("GTRE_Z","TRE_Z")){beta_0  <- beta_0_st + exp_u}
  
if(model_name %in% c("GTRE","TRE")){
if (isTRUE(is.numeric(start_val))){start_v <- start_val}else{
start_v <- if(is.na(beta_0_st)) {unname(c(lambda,sigma,sigma_r,sigma_h,beta_hat))} else{unname(c(lambda,sigma,sigma_r,sigma_h,beta_0,beta_hat))}}

out           <- matrix(0,nrow = 3,ncol = length(start_v))
rownames(out) <- c("par","st_err","t-val") 
if(model_name == "TRE"  & isTRUE(is.numeric(start_val))){colnames(out)<- c("lambda","sig","sig_r",x_vars_vec)}else{
colnames(out) <- c("lambda","sig","sig_r","sig_h",x_vars_vec)}

if(model_name == "TRE" & isFALSE(is.numeric(start_val))){out      <- out[,   -c(4)]  
                                                         start_v  <- start_v[-c(4)]}}  
}

results <- list(sfa_eps, sfa_alp, exp_eta, exp_u, sigma_v, sigma_u, sigma_r,
                sigma_h, beta_0, lambda, sigma, start_v, out,
                plm_gtre, beta_hat, alpha_hat, epsilon_hat, beta_0_st, beta_se,
                plm_tfe, plm_fd)  
  
names(results)  <- c("sfa_eps", "sfa_alp", "exp_eta", "exp_u", "sigma_v", "sigma_u", "sigma_r",
                     "sigma_h", "beta_0", "lambda", "sigma", "start_v", "out", 
                     "plm_gtre", "beta_hat", "alpha_hat", "epsilon_hat", "beta_0_st", "beta_se",
                     "plm_tfe", "plm_fd")
return(results)
}

start.tfe  <- function(formula_x, data, model_name, start_val, intercept, x_vars_vec, gamma, individual, N, y_var, n_x_vars){  
plm_tfe         <- plm(formula_x,  data,effect = "individual",model = "within" )
if (isTRUE(is.numeric(start_val))) {start_v <- start_val}else{ 
beta_hat         <- plm_tfe$coefficients[c(x_vars_vec)]
epsilon_hat      <- plm_tfe$residuals
beta_se          <- as.data.frame(summary(plm_tfe)[1])$coefficients.Std..Error[-c(1)]
sfa_eps          <- pcs_c(Y  = as.numeric(epsilon_hat))[[1]]$par
exp_u            <- sfa_eps[3]
sigma_v          <- sqrt(sfa_eps[2]^2/ (1 + sfa_eps[1]^2))    
sigma_u          <- sigma_v*sfa_eps[1]                      
lambda           <- sigma_u/sigma_v
sigma            <- sqrt(sigma_u^2 + sigma_v^2)
start_v          <- unname(c(lambda,sigma,beta_hat))
start_v[1]       <- if(gamma==TRUE) {sigma_u^2/sigma^2} else{start_v[1]} }
  
out<- matrix(0,nrow = 3,ncol = length(start_v))
rownames(out)  <- c("par","st_err","t-val") 
colnames(out)  <- c("lambda","sig",c(names(plm_tfe$coefficients))  ) 
colnames(out)[1] <- if(gamma==TRUE) {"gamma"} else{"lambda"}
  
indiv   <- noquote(as.vector(unique(data[,c(individual)])))
t       <- rep(0,N)
data_i  <- Y <- eps <- data_i_vars <- one_t <- I_t <- one_t1  <- I_t1  <- as.list(rep(0,N))
  
for (ii in 1:N) {
data_i[[ii]]   <-  data[which(data[,c(individual)]==indiv[ii]),]
t[ii]          <-  nrow(data_i[[ii]])
Y[[ii]]        <-  demean(as.numeric(data_i[[ii]][,y_var]))
data_i_vars[[ii]]   <- data.frame(data_i[[ii]][,c(x_vars_vec)] )
one_t[[ii]]    <- rep(1,t[[ii]])
I_t[[ii]]      <- diag(t[[ii]]) 
one_t1[[ii]]   <- rep(1,t[[ii]]-1)
I_t1[[ii]]     <- diag(t[[ii]]-1)  }
  
if(gamma==TRUE){upper <- c(1,rep(Inf,n_x_vars+1))} else{upper <- NA}  
  
results        <-  list(upper, I_t1, one_t1, I_t, one_t, data_i_vars, Y, t, data_i,
                     eps, indiv, out, start_v)
names(results) <- c("upper", "I_t1", "one_t1", "I_t", "one_t", "data_i_vars", "Y", "t", "data_i",
                           "eps", "indiv", "out", "start_v")
return(results)}

lower.start <- function(start_v, model_name, differ){
if(model_name == "TRE"){                                                 lower1 <- c(rep(.0000001,3) ,  start_v[-c(1:3)] - differ)}
if(model_name == "GTRE"){                                                lower1 <- c(rep(.0000001,4) ,  start_v[-c(1:4)] - differ)} 
if(model_name %in% c("NHN","NE","NR","NTN")){                            lower1 <- c(rep(.0000001,2),   start_v[-c(1:2)] - differ )}
if(model_name =="THT"){                                                  lower1 <- c(rep(.0000001,3),   start_v[-c(1:3)] - differ )}  
if(model_name %in% c("NG","NNAK")){                                      lower1 <- c(rep(.0000001,3),   start_v[-c(1:3)] - differ )} 
if(model_name %in% c("ZISF")){                                           lower1 <- c(start_v[1]-differ, rep(.0000001,2),   start_v[-c(1:3)] - differ)}  
if(model_name %in% c("ZISF_Z")){                                         lower1 <- c(rep(.0000001,2),  start_v[-c(1:2)] - differ)}
if(model_name %in% c("NHN_Z","NE_Z")){                                   lower1 <- c(rep(.0000001,1),  start_v[-c(1)]   - differ)}  
upper1 <-  c(start_v+differ)  

results <- list(upper1,lower1)
names(results) <- c("upper1", "lower1")
return(results)}

start.time  <- function(){
start_time <- Sys.time()
return(start_time)}

end.time    <- function(start_time){
end_time   <- Sys.time()
total_time <- end_time - start_time
return(total_time)}
