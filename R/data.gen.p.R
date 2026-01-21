data_gen_p <-
function(t, N, rand, sig_u, sig_v, sig_r, sig_h, cons, tau=0.5, mu=0, beta1, beta2){

if(!is.null(rand)){
set.seed(rand)}    
    
    n         <- N*t                                 ## Observations
    x1        <- log(runif(n,0,10))                  ## x1 variable
    x2        <- log(runif(n,1,50))                  ## x2 variable
    r         <- rep(rnorm(N,0,sig_r),each=t)        ## r_i individual effects
    u         <- abs(rnorm(n,0,sig_u))               ## u_it half-normal error term
    v         <- rnorm(n,0, sig_v)                   ## v_it error term 
    h         <- abs(rep(rnorm(N,0,sig_h),each=t))   ## h_i individual effects
    name      <- rep(1:N,each=t)                     ## name for individuals 
    year      <- rep(seq(1,t,1),N)                   ## general time frame
    cons      <- cons                                ## constant
    tau       <- tau                                 ## dependence paramter for tfe 
    x1_w      <- tau*r + sqrt(1-tau^2)* x1           ## x1 variable with dependence 
    x2_w      <- tau*r + sqrt(1-tau^2)* x2           ## x2 variable with dependence 
    
    
    ## make a first-difference data set from WangHo2010
    r_fd      <- runif(N,min = 0, max = 1)
    x_fd      <- rep(0,n)
    
    for (i in 1:N) {
      B          <- t*i
      A          <- B - (t - 1)
      x_fd[A:B]  <- rnorm(t, mean = r_fd[i], sd= 1)}
    
    z_fd       <- rnorm(n,0,1) 
    u_fd_star  <- abs(rep(rnorm(N, mean = mu, sd=sig_u),each=t))
    r_fd       <- rep(r_fd, each = t)
    u_fd       <- exp(mu*z_fd) * u_fd_star
    
    ## Output -  psfm 
    y_gtre    <- r - h + cons + beta1*x1    + beta2*x2   + v - u
    y_tre     <- r     + cons + beta1*x1    + beta2*x2   + v - u
    y_pcs     <-         cons + beta1*x1    + beta2*x2   + v - u
    y_tfe    <- r            + beta1*x1_w  + beta2*x2_w + v - u
    y_fd      <- r_fd         + beta1*x_fd               + v - u_fd
    
    y_gtre_nc <- r - h +        beta1*x1    + beta2*x2   + v - u
    y_tre_nc  <- r     +        beta1*x1    + beta2*x2   + v - u
    y_pcs_nc  <-                beta1*x1    + beta2*x2   + v - u
    
    c_gtre    <- r + h + cons + beta1*x1    + beta2*x2   + v + u
    c_tre     <- r     + cons + beta1*x1    + beta2*x2   + v + u
    c_pcs     <-         cons + beta1*x1    + beta2*x2   + v + u
    c_tfe    <- r            + beta1*x1_w  + beta2*x2_w + v + u
    
    c_gtre_nc <- r + h +        beta1*x1    + beta2*x2   + v + u
    c_tre_nc  <- r     +        beta1*x1    + beta2*x2   + v + u
    c_pcs_nc  <-                beta1*x1    + beta2*x2   + v + u
    
    ## GTRE-Z
    u_gtre         <- rep(0,n)
    z_gtre         <- runif(n,1,2)
    ## u_it half-normal error term with z var and cons
    for (i in 1:n) {
      u_gtre[i]      <- abs(rnorm(1,0, exp(0.9 + 0.6*z_gtre[i] )))   
    }
    
    y_gtre_z   <-         r  - h  +  cons + beta1*x1    + beta2*x2   + v - u_gtre
    y_tre_z   <-          r       +  cons + beta1*x1    + beta2*x2   + v - u_gtre
    
    data_trial <- as.data.frame(cbind(name,  year,  cons,  x1,  x1_w,  x2,  x2_w,  u,  v,  r,  y_gtre,  y_tre,  y_tfe,  h,  y_gtre_nc,  y_tre_nc,  y_pcs,  y_pcs_nc,  c_gtre,  c_tre, c_tfe,   c_gtre_nc,   c_tre_nc,   c_pcs,   c_pcs_nc,  y_fd,   x_fd,  u_fd_star,   z_fd,   r_fd,   u_fd,  u_gtre,  z_gtre,  y_gtre_z,   y_tre_z )) 
    colnames(data_trial) <- c(       "name","year","cons","x1","x1_w","x2","x2_w","u","v","r","y_gtre","y_tre","y_tfe","h","y_gtre_nc","y_tre_nc","y_pcs","y_pcs_nc","c_gtre","c_tre","c_tfe","c_gtre_nc", "c_tre_nc", "c_pcs", "c_pcs_nc","y_fd", "x_fd","u_fd_star", "z_fd", "r_fd", "u_fd","u_gtre","z_gtre","y_gtre_z", "y_tre_z")
    data_rand           <- pdata.frame(data_trial,    c("name","year"))
    
    return(data_rand) 
  }
