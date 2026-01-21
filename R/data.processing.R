# Package-level constants
.SFA_CONSTANTS <- list(
  CLIP_Z1_UPPER  =  37,
  CLIP_Z1_LOWER  = -37,
  CLIP_Z2_UPPER  =   8,
  CLIP_Z2_LOWER  = -37,
  HALTON_DISCARD =  1000,
  MIN_POSITIVE   = .Machine$double.eps,
  MAX_VALUE      = .Machine$double.xmax
  # HALTON_PRIMES = c(2, 3)
  )

## Second data cleaning for missing values  
data_proc2 <- function(data, data_x, fancy_vars, fancy_vars_z, data_z, y_var, 
                       x_vars_vec, halton_num, individual, N, model_name, rand.gtre){

data           <- data[rownames(data_x),]
  
if (isTRUE(length(fancy_vars)>0)) { data_inter   <- cbind(  data, data_x[,fancy_vars]  )
colnames(data_inter) <- c(colnames(data),fancy_vars)
data  <- data_inter} else {}
  
if (isTRUE(length(fancy_vars_z)>0)) { data_inter   <- cbind(  data, data_z[,fancy_vars_z]  )
colnames(data_inter) <- c(colnames(data),fancy_vars_z)
data  <- data_inter} else {}  
  
Y             <-  as.matrix(subset(data,select = y_var))
data_i_vars   <-  as.matrix(data.frame(subset(data,select = x_vars_vec)))

if(model_name %in% c("GTRE","TRE") ){
if (isTRUE(is.numeric(halton_num))) {R <- halton_num}else{
R       <- ceiling(sqrt(nrow(data)))+100 }         ## Integral reps  

  R_H     <- randtoolbox::halton(R+1000,2,start = 1,normal = FALSE)[-c(1:1000),c(1:2)]   
  R_H     <- cbind( qnorm(R_H[,1]) , sqrt(2)* pracma::erfinv(R_H[,2]) )                ## using inverse error function for R_H2
  
  if(!is.null(rand.gtre)){
     set.seed(rand.gtre)}   
  
  mat <- matrix(0,nrow=R, ncol=9999)
  for(v in 1:9999){mat[,v] <-  sample(R_H[,1])}
  
  cor <- matrix(0,9999,1)
  for(v in 1:9999){cor[v] <-  abs(cor(mat[,v],R_H[,2]))}
  
  R_H <- cbind(mat[,which(cor==min(cor))] , R_H[,2])
  rm(cor,v,mat)

# print(paste( "Primes 2 and 3 are in use, with 1,000 discards.  Correlation between R and H draws is:", round(cor(R_H)[1,2],10), sep = "" )) 

indiv   <- noquote(as.vector(unique(data[,c(individual)])))
t       <- rep(0, N)
data_i  <- Y <- eps <- data_i_vars <- R_h1 <- R_h2 <- as.list(rep(0,N))
  
for (ii in 1:N){
data_i[[ii]]  <-  data[which(data[,c(individual)]==indiv[ii]),]
t[ii]         <-  nrow(data_i[[ii]])
R_h1[[ii]]    <-  t(matrix(rep(R_H[,1],t[[ii]]),R,t[[ii]]))
R_h2[[ii]]    <-  abs(t(matrix(rep(R_H[,2],t[[ii]]),R,t[[ii]])))
Y[[ii]]       <-  matrix(rep(data_i[[ii]][,y_var],R),t[[ii]],R)
data_i_vars[[ii]]   <- data.frame(data_i[[ii]][,c(x_vars_vec)] )}  

result_gtre_tre         <- list(R,   R_H,   indiv,   t,   data_i,   eps,   R_h1,   R_h2,  data_x)
names(result_gtre_tre)  <-  c( "R", "R_H", "indiv", "t", "data_i", "eps", "R_h1", "R_h2","data_x")
}  

result         <- list(data, Y,  data_i_vars)
names(result)  <- c("data", "Y","data_i_vars")

if (model_name %in% c("GTRE", "TRE")) {result <- c(result, result_gtre_tre)}

return(result)

}

## Inital Data cleaning and Processing 
data_proc <- function(formula, data, model_name, individual=NULL, inefdec){
data_z     <- NULL
data_orig  <- data
form_parts <- base::strsplit(as.character(formula), "|", fixed = TRUE)

if(isFALSE(unique(grepl( "|", deparse(formula), fixed = TRUE))) ){formula_x  <- formula
y_var     <- gsub(" ", "", noquote(as.character( unlist( strsplit( deparse(formula_x), "~", fixed = TRUE)[[1]][1]))))}else{
formula_x <- paste(form_parts[[2]], "~",form_parts[[3]][1], sep = "")
y_var     <- gsub(" ", "", noquote(as.character( unlist(strsplit( formula_x, "~", fixed = TRUE))[[1]])))}

if(length(unlist(form_parts))>3){
formula_z <- paste(form_parts[[2]], "~",form_parts[[3]][2], sep = "")
z_vars    <- as.character(gsub(" ", "", noquote(as.character(unlist(form_parts)[[4]])))) 
z_vars    <- noquote(gsub("+", " ", z_vars, fixed=TRUE))
z_vars    <- unlist(strsplit(z_vars, " "))
data_z    <- data_conform(formula = formula_z, data = data)
if(model_name=="NHN"){ model_name <- "NHN_Z"}
if(model_name=="NE"){  model_name <- "NE_Z"}
if(model_name=="GTRE"){model_name <- "GTRE_Z"}
if(model_name=="TRE"){ model_name <- "TRE_Z"}
if(model_name %in% c("THT","NTN","CHC","NU")){return(c("Currently building this functionality"))}}#else{z_vars <- NULL}

data_x    <- data_conform(formula = formula_x, data = data)

intercept      <- if(isTRUE(grepl(-1, gsub("[[:space:]]", "",as.character(formula_x))))) {0} else{1}
inefdec_n      <- if(isTRUE(inefdec) ) {1} else{-1}
inefdec_TF     <- if(isTRUE(inefdec) ) {TRUE} else{FALSE}
x_vars_vec     <- if(model_name %in% c("TFE","FD") & intercept==1){colnames(data_x)[-c(1)]}else {colnames(data_x)}
n_x_vars       <- length(x_vars_vec)
x_vars         <- x_vars_vec
x_x_vec        <- rep(0,length= n_x_vars)
fancy_vars     <- setdiff(colnames(data_x),colnames(data))
fancy_vars_z   <- NULL 

N <- NULL 
if(isTRUE(individual==NULL)){}else{N  <- length(unique(data[,c(individual)])) }

n_z_vars       <- NA  ## maybe make condition here making this default to null
formula_z <- z_vars <- intercept_z <- z_vars_vec <- z_z_vec <- n_z_vars <- NULL

if(length(unlist(form_parts))>3){    
intercept_z    <- if(isTRUE(grepl(-1, gsub("[[:space:]]", "",as.character(formula_z))))) {0} else{1}
z_vars_vec     <- if(model_name %in% c("TFE","FD") & intercept_z==1){colnames(data_z)[-c(1)]}else {colnames(data_z)}
n_z_vars       <- length(z_vars_vec)
z_vars         <- z_vars_vec
z_z_vec        <- rep(0,length= n_z_vars)
fancy_vars_z   <- setdiff(colnames(data_z),colnames(data))
}

result  <-  list(data_orig, form_parts, formula_x, y_var, model_name, data_x, intercept, inefdec_n, inefdec_TF,
                                       x_vars_vec, n_x_vars, x_vars, x_x_vec, fancy_vars,fancy_vars_z,n_z_vars,N,
                 formula_z, z_vars, data_z, intercept_z, z_vars_vec, z_z_vec, n_z_vars, data_z)

names(result)<- c("data_orig", "form_parts", "formula_x", "y_var", "model_name","data_x", "intercept", "inefdec_n",
                  "inefdec_TF","x_vars_vec", "n_x_vars", "x_vars", "x_x_vec", "fancy_vars","fancy_vars_z","n_z_vars","N",
                  "formula_z", "z_vars", "data_z", "intercept_z", "z_vars_vec", "z_z_vec", "n_z_vars","data_z")
return(result)
}

## lm code for data   
data_conform <- function (formula, data,na.action, method = "qr", 
                          model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
                          contrasts = NULL, offset, ...){
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (method == "model.frame")return(mf) else if (method != "qr") 
    warning(gettextf("method = '%s' is not supported. Using 'qr'", method), domain = NA)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf, contrasts)
  return(x)}

## Basic PCS code for regression on intercept
pcs_c  <- function(Y, inefdec=TRUE, Method= formals(psfm)$Method){
sigma_u        <- runif(n=1, min = 0.04, max = 1) ## Random starting values
sigma_v        <- runif(n=1, min = 0.04, max = 1) 
beta_0         <- runif(n=1, min = 0.04, max = 1) 

lambda         <- sigma_u/sigma_v
sigma          <- sqrt(sigma_u^2 + sigma_v^2)
  
start_v        <- unname(c(lambda,sigma,beta_0))
inefdec_n      <- if(isTRUE(inefdec) ){1}else{-1}
  
fn = function(x){
eps     <- inefdec_n*(Y - x[3])
like  <- sum( log( (2/x[2])  * 
             pmax(dnorm(eps/x[2]),eps*0+.Machine$double.eps)*  
             pmax(pnorm(-eps*x[1]/x[2]),eps*0+.Machine$double.eps)   ))
return(-like[is.finite(like)])}  
opt <- optim(par = start_v, 
             fn = fn,
             lower = c(rep(0.01,2),-Inf),
             control = list(maxit = 300, REPORT = 1, trace = 0, pgtol=0),
             hessian = FALSE, 
             method = Method)
return(list(opt, Y))}


summary.sfareg <- function(object, ...){  
  # Perform specific summary calculations and formatting
  cat("--- SFA Regression Model Summary ---\n")
  cat("Formula:", deparse(object$formula), "\n")
  cat("Total time:", object$total_time, "\n")
  cat("Model Output:\n")
  print(object$out)
  
  ## return original object 
  invisible(object) }


print.sfareg <- function(x, ...) {
  # Perform specific summary calculations and formatting
  cat("--- SFA Regression Model Summary ---\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Total time:", x$total_time, "\n")
  cat("Model Output:\n")
  print(x$out)}

