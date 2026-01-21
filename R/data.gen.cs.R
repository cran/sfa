data_gen_cs <- function(N, rand, sig_u, sig_v, cons, beta1, beta2, a, mu){
if(!is.null(rand)){
set.seed(rand)}    

n        <- N                                   ## Observations
x1       <- log(runif(n,2,10))                  ## x1 variable
x2       <- log(runif(n,3,50))                  ## x2 variable
    
## Normal half normal
u        <- abs(rnorm(n,0,sig_u))               ## u_it half-normal error term
v        <- rnorm(    n,0,sig_v)                
    
## Normal truncated normal
u_tn     <- rtruncnorm(n=n, mean = mu, sd=sig_u, a=0)

## NHN with Z
uz         <- rep(0,n)
z          <- runif(n,1,2)
    
for (i in 1:n) {
uz[i]      <- abs(rnorm(1,0, exp(0.9 + 0.6*z[i] )  )) ## u_it half-normal error term with z var and cons  
}

## Normal exponetial
phi      <- 1/sig_u
u_e      <- rexp(n,rate=phi)                    ## u_it exponential
    
## T half T
u_t      <- sig_u*abs(rt(n,a))                  ## u_it half-T error term
v_t      <- sig_v*    rt(n,a)                   
    
## Cauchy half Cauchy
u_c      <- abs(rcauchy(n,0,sig_u))             ## u_it half-cauchy error term
v_c      <- rcauchy(    n,0,sig_v)             
    
## Normal half uniform
u_u      <- runif(n, min=0, max=sig_u)          ## u_it half-uniform error term
name     <- rep(1:N,each=t)                     ## name for individuals 
cons     <- cons                                ## constant
    
## Output - sfm 
y_pcs    <-         cons + beta1*x1    + beta2*x2   + v   - u
y_pcs_z  <-         cons + beta1*x1    + beta2*x2   + v   - uz
y_pcs_t  <-         cons + beta1*x1    + beta2*x2   + v_t - u_t
y_pcs_tn <-         cons + beta1*x1    + beta2*x2   + v   - u_tn
y_pcs_e  <-         cons + beta1*x1    + beta2*x2   + v   - u_e
y_pcs_c  <-         cons + beta1*x1    + beta2*x2   + v_c - u_c
y_pcs_u  <-         cons + beta1*x1    + beta2*x2   + v   - u_u
    
## Normal + Cauchy - half-normal
y_pcs_w  <-         cons + beta1*x1    + beta2*x2   + v +v_c   - u
    
data_trial <- as.data.frame(cbind(name,  cons,  x1,  x2,  u, uz,   v,  u_t,  v_t,  u_c,  v_c,  u_e , u_u,   y_pcs,  y_pcs_t, y_pcs_e,   y_pcs_c,  y_pcs_u,   y_pcs_z,   y_pcs_w,   y_pcs_tn,   z))                     
colnames(data_trial) <- c(       "name","cons","x1","x2","u","uz","v","u_t","v_t","u_c","v_c","u_e","u_u" ,"y_pcs","y_pcs_t","y_pcs_e","y_pcs_c","y_pcs_u", "y_pcs_z", "y_pcs_w" ,"y_pcs_tn", "z")
    
data_rand           <- data.frame(data_trial)
data_rand$con       <- 1
return(data_rand) 
}