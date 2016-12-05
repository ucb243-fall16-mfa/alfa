if("ggplot2" %in% rownames(installed.packages()) == FALSE) {
  install.packages("ggplot2")}
if("grid" %in% rownames(installed.packages()) == FALSE) {
  install.packages("grid")}
if("gridExtra" %in% rownames(installed.packages()) == FALSE) {
  install.packages("gridExtra")}

# ---------------------------------
# General check functions
# ---------------------------------

check_data <- function(data) {
  if (class(data)!="data.frame" & class(data)!="matrix") {
    stop("Data should be data frame or matrix")
  }
  if (sum(!sapply(data,is.numeric))!=0) {
    stop("Data is not numeric")
  }
  TRUE
}

#### more to test
check_sets <- function(sets) {
  if(is.list(sets)==FALSE){
    stop("sets should be a list")
  }
  if (!sum(sapply(sets,is.numeric))!= 0){
    TRUE
  }else if (!sum(sapply(sets,is.character))!= 0){
    TRUE
  }else{
    stop("sets should be all numeric or all character")
  }
}


check_ncomps <- function(data, ncomps) {
  if (!is.null(ncomps) & !is.numeric(ncomps)) {
    stop("ncomps should be a numeric number!")
  }
  if (!is.null(ncomps) & length(ncomps) != 1) {
    stop("The length of ncomps should be 1!")
  }
  if (!is.null(ncomps)) {
    if (!ncomps %in% seq(from = 1, to = nrow(data), by = 1)) {
      stop(paste("The value of ncomps should be between 1 and", nrow(data)))
    }
  }
  TRUE
}

#check numeric
check_numeric <- function(x){
  if(sum(!sapply(x,is.numeric))!=0){
    stop("Argument is not numeric")
  }
  TRUE
}



# ------------------------------------

#' @title Construct mfa
#' @description Creates an object of class \code{"mfa"}
#' @param data dataset on which MFA is conducted
#' @param sets list of vectors indicating the sets of variables (i.e. the blocks)
#' @param ncomps integer indicating how many number of components (i.e. factors) are to be extracted
#' @param center either a logical value or a numeric vector of length equal to the number of active variables in the analysis
#' @param scale either a logical value or a numeric vector of length equal to the number of active variables in the analysis
#' @return an object of class mfa
#' @export
#' @examples
#' url="https://raw.githubusercontent.com/ucb-stat243/stat243-fall-2016/master/problem-sets/final-project/data/wines.csv"
#' data=read.csv(file=url, header = TRUE, sep = ",")
#' datanum <- data[,-1]
#'
#' # a list of numeric vectors with the position of the active variables in the data table.
#' sets <- list(1:6, 7:12, 13:18, 19:23, 24:29, 30:34, 35:38, 39:44, 45:49, 50:53)
#' # create object 'wine' of class mfa
#' wine <- mfa_gen(data = datanum, sets = sets, ncomps = 2)

mfa_gen <- function(data, sets, ncomps = NULL, center = TRUE, scale = TRUE) {
  check_data(data)
  check_sets(sets)
  check_ncomps(data, ncomps)
  #set ncomps
  ncomps <- ifelse(is.null(ncomps), nrow(data), ncomps)


  # store the weights
  alpha <- NULL

  # data.pro stores the normalized grand table
  data.pro <- matrix(nrow=nrow(data), ncol=tail(sets[[length(sets)]],n=1))

  # conduct SVD on each sub table
  for (i in 1:length(sets)){
    #decompose the grand table
    dat <- data[,sets[[i]]]
    #normalize subtable
    if (all(scale==F)){
      dat <- scale(dat,center=center,scale=F)
    }else if(all(scale==T)){
      dat <- scale(dat,center=center,scale=T)/sqrt(nrow(data)-1)
    }else{
      dat <- scale(dat,center=center, scale=scale[sets[[i]]])/sqrt(nrow(data)-1)
    }
    
    data.pro[,sets[[i]]] <- dat

    # SVD
    s <- svd(dat)
    alpha[i] <- 1/(s$d[1]^2)
  }

  # diagonal matrix
  A <- diag(rep(alpha,lapply(sets,length)))


  #compute m, p, q for gsvd
  M <- diag(1/nrow(data),nrow=nrow(data))
  gsvd <- svd(sqrt(M)%*%data.pro%*%sqrt(A))
  P <- solve(sqrt(M))%*%gsvd$u
  Q <- solve(sqrt(t(A)))%*%gsvd$v
  D <- diag(gsvd$d)


  # Full Factor Scores
  Fscore <- P %*% D
  # requested Factor Scores
  Factorscore = P%*%D[,1:ncomps]


  # partial factor scores
  # create a list of matrix to store the PFS for each experts
  partial.factor <- vector("list", length(sets))
  for (i in 1:length(sets)){
    #decompose the processed data table
    dat <- data.pro[,sets[[i]]]
    partial.factor[[i]] <- length(sets)*alpha[i]*dat%*%Q[sets[[i]],1:ncomps]
  }


  # vector of eigenvalues
  eigenvalues <- eigenvalues <- (diag(D)[1:ncomps])^2

  # create the object mfa
  object <- list("eigenvalues" = eigenvalues, "FactorScore" = Factorscore,
                 "PartialFactorScores" = partial.factor, "Loadings" = Q[,1:ncomps],
                 "MatrixA" = A, "sets" = sets)
  class(object) <- "mfa"
  return(object)
}


# ===============================
# print method
# ===============================

# ---------------------------------
# Auxiliary check function for print method
# ---------------------------------
check_subtable <- function(mfa,subtable) {
  if (!is.numeric(subtable) & !is.null(subtable)) {
    stop("subtable should be either a numeric value or NULL")
  }
  if(!is.null(subtable)) {
    if(!subtable %in% 1:length(mfa$sets)){
      stop(paste("subtable should be between 1 and", length(mfa$sets)))
    }
  }
  TRUE
}

# ---------------------------------
#' @title Print method for object of class \code{"mfa"}
#' @description Print out the eigenvalues and partial factor score of the specified subtable
#' @param object mfa object
#' @param subtable if subtable = NULL then no partial factor score will be printed; if subtable is a numerical value, then the partical factor score for the subtable of this index is printed
#' @return Printed eigenvalues and partical factor score
#' @export
#' @examples
#' # print the eigenvalues
#' print(wine)
#'
#' # print the eigenvalues and the partical factor scores for the 6th subtable
#' print(wine, subtable = 6)

print.mfa <- function(object, subtable = NULL) {

  check_subtable(object, subtable)
  # print eignvalues
  cat('Eigenvalues (listed in the order of components) \n')
  e = matrix(round(object$eigenvalues, 3), nrow = 1, ncol = length(object$eigenvalues))
  # change the row names and column names of the matrix
  rownames(e) = ""
  colnames(e) = 1:length(object$eigenvalues)
  print(e)

  cat("\n")

  # print partial factor scores if subtable index is indicated
  if (class(subtable) == "numeric") {
    cat(paste('Partial Factor Score Matrix for Assessor No.', subtable, '\n', sep = ""))

    # change the row names and column names of the matrix
    rownames(object$`PartialFactorScores`[[subtable]]) = 1:nrow(object$FactorScore)
    colnames(object$`PartialFactorScores`[[subtable]]) = 1:ncol(object$Loadings)
    # print the matrix
    print(object$`PartialFactorScores`[[subtable]])
  }

  # hide the default description of object
  invisible(object)
}

# ===============================
# summaries of eigenvalues
# ===============================

#' @title Eigenvalue summary for object of class \code{"mfa"}
#' @description Creates a table that summarizes information about the obtained eigenvalues
#' @param x mfa object
#' @return a table with singular values, the eigenvalues, cumulative, percentage of intertia, cumulative percentage of inertia, for all the extracted components.
#' @export
#' @examples
#' eigenval(wine)

eigenval <- function(x) UseMethod("eigenval")

#' @rdname eigenval
#' @export
eigenval.mfa <- function(x){
  eigenvalue <- round(as.numeric(x$eigenvalues),3)
  singular <- round(sqrt(eigenvalue),3)
  cumulative <- cumsum(eigenvalue)
  intertia <- round(eigenvalue/sum(eigenvalue)*100,0)
  cumulative.intertia <- cumsum(intertia)
  eigenvalue.table <- t(data.frame(singular,eigenvalue,cumulative, intertia, cumulative.intertia))
  colnames(eigenvalue.table) <- 1:length(eigenvalue)
  return(eigenvalue.table)
}

# ===============================
# contributions
# ===============================

#' @title Contribution of observations
#' @description Calculates the contributions of the observations to the extracted components (dimensions)
#' @param x mfa object
#' @return a matrix with its [i,j] element be the contribution of the ith observation to the jth component
#' @export
#' @examples
#' contri_obs(wine)

# Define the generic function
contri_obs <- function(x) UseMethod("contri_obs")

#' @rdname contri_obs
#' @export
contri_obs.mfa <- function(x) {
  # Get the eigenvalues
  lambda <- x$eigenvalues[1:ncol(x$Loadings)]
  # Inverse the eigenvalues
  inv_lambda <- lambda^(-1)
  # Matrix of mass
  M <- diag(1/nrow(x$FactorScore),nrow=nrow(x$FactorScore))
  # Square the factor scores
  F_squared <- (x$FactorScore)^2
  # Contribution of an observation to a dimension
  ctr_obs <- M %*% F_squared * inv_lambda

  return(ctr_obs)
}

#' @title Contribution of variables
#' @description Calculates the contributions of the variables to the extracted components (dimensions)
#' @param x mfa object
#' @return a matrix with its [i,j] element be the contribution of the ith variable to the jth component
#' @export
#' @examples
#' contri_var(wine)

# Define the generic function
contri_var <- function(x) UseMethod("contri_var")

#' @rdname contri_var
#' @export
contri_var.mfa <- function(x) {
  # Matrix A
  A <- x$MatrixA
  # Square the loadings
  Q_squared <- (x$Loadings)^2
  # Contribution of a variable to a dimension
  ctr_var <- A %*% Q_squared

  return (ctr_var)
}

#' @title Contribution of subtables
#' @description Calculates the contributions of the subtables to the extracted components (dimensions)
#' @param x mfa object
#' @return a matrix with its [i,j] element be the contribution of the ith subtable to the jth component
#' @export
#' @examples
#' contri_table(wine)

# Define the generic function
contri_table <- function(x) UseMethod("contri_table")

#' @rdname contri_table
#' @export
contri_table.mfa <- function(x) {
  # Matrix A
  A <- x$MatrixA
  # Square the loadings
  Q_squared <- (x$Loadings)^2
  # Contribution of a variable to a dimension
  ctr_var <- A %*% Q_squared
  # Get the sets
  set <- x$sets
  ctr_table <- matrix(nrow = length(set), ncol = ncol(ctr_var))
  for (i in 1:length(set)) {
    ctr_table[i,] <- colSums(ctr_var[set[[i]],])
  }
  return(ctr_table)
}

# ===============================
# Rv coefficient
# ===============================

#' @title Rv coefficient of two tables
#' @description Calculates the Rv coefficients between 2 tables
#' @param table1 data table 
#' @param table2 data table
#' @return the Rv coefficient between table1 and table2
#' @export
#' @examples
#' t1 <- matrix(c(2, 4, 3, 1, 5, 7))
#' t2 <- matrix(c(1, 8, 9, 10, 4, 6))
#' RV(table1 = t1, table2 = t2)
RV <- function(table1,table2){
  check_numeric(table1)
  check_numeric(table2)
  table1 <- as.matrix(table1)
  table2 <- as.matrix(table2)
  numerator <- sum(diag((table1 %*% t(table1)) %*% (table2 %*% t(table2))))
  denominator <- sqrt((sum(diag((table1 %*% t(table1)) %*% (table1 %*% t(table1))))) * (sum(diag((table2 %*% t(table2)) %*% (table2 %*% t(table2))))))
  return(numerator/denominator)
}


#' @title Rv coefficient of a dataset
#' @description Calculates the Rv coefficients between the specified subtables in a grand table
#' @param dataset the grand table from which subtables are extracted
#' @param sets list of numeric vectors that extract subtables from the grand table. The ith vector of the list represents the ith subtable extracted
#' @return a matrix with its [i,j] element be the Rv coefficient between the ith subtable and the jth subtable
#' @export
#' @examples
#' RV_table(datanum, sets = list(sets[[1]], sets[[2]]))

RV_table <- function(dataset, sets){

  check_data(dataset)
  check_sets(sets)

  k <- length(sets)
  table <- matrix(nrow= k, ncol= k)
  for (i in 1:k){
    for (j in 1:k){
      table[i,j] <- RV(dataset[,sets[[i]]], dataset[,sets[[j]]])
    }
  }
  return(table)
}

# ===============================
# Lg coefficient
# ===============================

# Auxiliary function to calculate Lg coefficient for 2 tables
Lg <- function(table1,table2,alpha1,alpha2){

  check_numeric(table1)
  check_numeric(table2)
  check_numeric(alpha1)
  check_numeric(alpha2)

  table1 <- as.matrix(table1)
  table2 <- as.matrix(table2)
  s <- sum(diag((table1 %*% t(table1)) %*% (table2 %*% t(table2))))
  return(s*alpha1*alpha2)
}

#' @title Lg coefficient
#' @description Calculates the Lg coefficients between the specified subtables in a grand table
#' @param dataset the grand table from which subtables are extracted
#' @param sets list of numeric vectors that extract subtables from the grand table. The ith vector of the list represents the ith subtable extracted
#' @return a matrix with its [i,j] element be the Lg coefficient between the ith subtable and the jth subtable
#' @export
#' @examples
#' Lg_table(datanum, sets = list(sets[[1]], sets[[2]]))

Lg_table <- function(dataset, sets){

  check_data(dataset)
  check_sets(sets)

  k <- length(sets)
  table <- matrix(nrow= k, ncol= k)
  for (i in 1:k){
    for (j in 1:k){
      gam1 <- svd(dataset[,sets[[i]]])$d[1]
      gam2 <- svd(dataset[,sets[[j]]])$d[1]
      a1 <- 1/(gam1^2)
      a2 <- 1/(gam2^2)
      table[i,j] <- Lg(dataset[,sets[[i]]], dataset[,sets[[j]]], alpha1 = a1, alpha2 = a2)
    }
  }
  return(table)
}

# ===============================
# Bootstrap
# ===============================

# Auxiliary function for bootstrapping one a time
bootstrapprep=function(dataset,userset,ncomp){
  bootPfactorscore=mfa_gen(dataset,sets=userset,ncomps =ncomp)$PartialFactorScores
  XB=sample(1:length(userset),length(userset),replace=TRUE)
  add=matrix(0,nrow = nrow(dataset),ncol=ncomp)
  for (i in 1:length(XB)){
    Factor=bootPfactorscore[[XB[i]]]
    add=add+Factor 
  }
  Fboot=1/length(XB)*add
  return (Fboot)
}


#' @title Bootstrap for object of class \code{"mfa"}
#' @description Perform bootstrap in order to get the matrix of bootstrap ratio. The ratios are used to find the observations that reliably contributes to a given component
#' @param L size of the bootstrap sample
#' @param userset list of vectors to partition the grand table into subtables
#' @param dataset the grand table
#' @param ncomps integer indicating how many number of components (i.e. factors) are to be extracted
#' @return matrix of bootstrap ratios
#' @export
#' @examples
#' sets <- list(1:6, 7:12, 13:18, 19:23, 24:29, 30:34, 35:38, 39:44, 45:49, 50:53)
#' bootstrap(1000,sets,datanum,2)

# bootstrap by L times
bootstrap=function(L,dataset,userset,ncomps){
  Fl=matrix(0,nrow = nrow(dataset),ncol=ncomps)
  Fstar=list()
  for (k in 1:L){
    Fstar[[k]]=bootstrapprep(dataset,userset,ncomps)
    Fl=Fl+Fstar[[k]]
  }
  #mean
  L=1000
  F.starbar=1/L*Fl
  #std 
  F.starstdsum=matrix(0,nrow = nrow(datanum),ncol=2)
  for (k in 1:L){
    F.starstdsum=F.starstdsum+(F.starbar-Fstar[[k]])^2
  }
  F.starstd=sqrt(1/L*F.starstdsum)
  #Ratio
  T.star=F.starbar * F.starstd^(-1)
  return (T.star)
}



# ---------------------------------
# Auxiliary check function for plot method
# ---------------------------------
check_dim <- function(mfa, dim1, dim2) {
  if (!is.numeric(dim1) | !is.numeric(dim2)) {
    stop("dim1 and dim2 should be numeric!")
  }
  if (length(dim1) != 1 | length(dim2) != 1) {
    stop("The length of dim1 and dim2 should be 1")
  }
  if (!dim1 %in% seq(from = 1, to = ncol(mfa$FactorScore), by = 1) | !dim2 %in% seq(from = 1, to = nrow(mfa$FactorScore), by = 1)) {
    stop(paste("dim1 and dim2 should be between 1 and", ncol(mfa$FactorScore)))
  }
  TRUE
}

check_type <- function(type) {
  if (!is.numeric(type)) {
    stop("type should be a numeric value")
  }

  if(length(type) != 1) {
    stop("The length of type should be 1!")
  }

  if (!type %in% c(1,2,3,4)) {
    stop("type should be numeric value 1, 2, 3, or 4!")
  }
  TRUE
}

check_text <- function(mfa, text) {
  if (!is.null(text) & length(text) != nrow(mfa$Loadings)) {
    stop(paste("The length of text should be", nrow(mfa$Loadings)))
  }
  TRUE
}

check_cat <- function(mfa, cat) {
  if (!is.null(cat)) {
    if(sum(sapply(cat, length)) != nrow(mfa$FactorScore)) {
      stop(paste("The total number of elements in cat should be", nrow(mfa$FactorScore)))
    }
  } 
  TRUE
}

check_var <- function(mfa, var) {
  if (!is.null(var)) {
    if (!is.numeric(var)) {
      stop ("var should be a vector containing numeric values")
    }
    
    if (sum(!var %in% 1:length(mfa$sets)) != 0) {
      stop(paste("var should contain integers between 1 and", length(mfa$sets)))
    }
  }
  TRUE
}


check_scale <- function(scale_x, scale_y) {
  if (!is.numeric(scale_x) | !is.numeric(scale_y)) {
    stop("scale_x and scale_y should be numeric!")
  }
  if (length(scale_x) != 1 | length(scale_y) != 1) {
    stop("The length of scale_x and scale_y should be 1")
  }
  TRUE
}

# ---------------------------------
# Auxiliary functions for plot method
# ---------------------------------
# only plot the partial factor score of one applicant

plot_fs <- function(mfa, dim1, dim2, cat = NULL) {
  fs <- as.data.frame(mfa$FactorScore)
  x <- fs[,dim1]
  y <- fs[,dim2]
  xlab <- paste0("comp", dim1)
  ylab <- paste0("comp", dim2)
  
  
  if(!is.null(cat)) {
    category <- NULL
    for (i in 1:length(cat)) {
      if (is.null(names(cat))) {
        category[cat[[i]]] <- rep(i, length(cat[i]))
      }
      else {
        category[cat[[i]]] <- names(cat)[i]
      }
    }
    
    
    fs$Category <- factor(category)
    
    
    # Plot factor scores
    ggplot() +
      geom_point(data = fs, aes(x = x, y = y, col = Category, group = Category)) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      labs(title = "Factor Scores", x = xlab, y = ylab)
  }
  
  else {
    ggplot() +
      geom_point(data = fs, aes(x = x, y = y)) +
      geom_text(data = fs, aes(x = x, y = y, label = rownames(fs)),hjust = 0, nudge_x = 0.05, size=3) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      labs(title = "Factor Scores", x = xlab, y = ylab)
  }
}


plot_pfs <- function(mfa, dim1, dim2, index, cat = NULL) {
  pfs <- as.data.frame(mfa$PartialFactorScores[[index]])
  x <- pfs[,dim1]
  y <- pfs[,dim2]
  xlab <- paste0("comp", dim1)
  ylab <- paste0("comp", dim2)
  title <- paste("Partial Factor Score for Subtable", index)
  
  if(!is.null(cat)) {
    indicator <- NULL
    for (i in 1:length(cat)) {
      if (is.null(names(cat))) {
        indicator[cat[[i]]] <- rep(i, length(cat[i]))
      }
      else {
        indicator[cat[[i]]] <- names(cat)[i]
      }
    }
    pfs$category <- factor(indicator)
    ggplot() +
      geom_point(data = pfs, aes(x = x, y = y, col = category, group = category)) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      theme(text = element_text(size=5))+
      labs(title = title, x = xlab, y = ylab) +
      scale_y_continuous(limits=c(min(y)-sd(y),
                                  max(y)+sd(y)))+
      scale_x_continuous(limits=c(min(x)-sd(x),
                                  max(x)+sd(x)))
  }
  
  else {
    ggplot() +
      geom_point(data = pfs, aes(x = x, y = y),size = 0.5) + 
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_text(data = pfs, aes(x = x, y = y, label = rownames(pfs)), hjust = 0, nudge_x = 0.05,size=2) +
      theme(text = element_text(size=5))+
      labs(title = title, x = xlab, y = ylab) +
      scale_y_continuous(limits=c(min(y)-sd(y),
                                  max(y)+sd(y)))+
      scale_x_continuous(limits=c(min(x)-sd(x),
                                  max(x)+sd(x)))
  }
}


# only plot the loadings of one applicant
plot_loadings <- function(mfa, dim1, dim2, index, scale_x, scale_y, text = NULL) {
  sets <- mfa$sets
  loadings <- as.data.frame(mfa$Loadings[sets[[index]],])
  xlab <- paste0("loading", dim1)
  ylab <- paste0("loading", dim2)
  title <- paste("Loadings for Subtable", index)
  
  if (!is.null(text)) {
    row.names(loadings) <- text[sets[[index]]]
  }
  loadings[,dim1] <- loadings[,dim1]*scale_x
  loadings[,dim2]<- loadings[,dim2]*scale_y
  ggplot(data = loadings, aes(x=loadings[,dim1], y=loadings[,dim2])) +
    geom_point(size=0.5) +
    geom_text(aes(label = rownames(loadings)), hjust = 0, nudge_x = 0.05,size=3) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    theme(text = element_text(size=5))+
    labs(title = title, x = xlab, y = ylab)+
    scale_y_continuous(limits=c(min(loadings[,dim2])-sd(loadings[,dim2]),
                                max(loadings[,dim2])+sd(loadings[,dim2])))+
    scale_x_continuous(limits=c(min(loadings[,dim1])-sd(loadings[,dim1]),
                                max(loadings[,dim1])+sd(loadings[,dim1])))
}


# plot the partial factor score and loadings of one applicant
plot_pfsl <- function(mfa, dim1, dim2, index, scale_x, scale_y, text = NULL, cat = NULL) {
  pfs <- as.data.frame(mfa$PartialFactorScores[[index]])
  x <- pfs[,dim1]
  y <- pfs[,dim2]
  xlab <- paste0("comp", dim1)
  ylab <- paste0("comp", dim2)
  title <- paste("Partial Factor Score and Loadings for Subtable", index)
  
  
  # loadings
  sets <- mfa$sets
  loadings <- as.data.frame(mfa$Loadings[sets[[index]],])
  loadings[,dim1] <- loadings[,dim1]*scale_x
  loadings[,dim2]<- loadings[,dim2]*scale_y
  if (!is.null(text)) {
    row.names(loadings) <- text[sets[[index]]]
  }
  
  
  # pfs
  if(!is.null(cat)) {
    indicator <- NULL
    for (i in 1:length(cat)) {
      if (is.null(names(cat))) {
        indicator[cat[[i]]] <- rep(i, length(cat[i]))
      }
      else {
        indicator[cat[[i]]] <- names(cat)[i]
      }
    }
    pfs$category <- factor(indicator)
    ggplot() +
      geom_point(data = pfs, aes(x = x, y = y, col = category, group = category),size=0.5) +
      geom_point(data = loadings, aes(x=loadings[,dim1], y=loadings[,dim2]),size=0.5) +
      geom_text(data = loadings, aes(x=loadings[,dim1], y=loadings[,dim2], label = rownames(loadings)), hjust = 0, nudge_x = 0.05, size=2) +
      geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
      theme(text = element_text(size=5))+
      labs(title = title, x = xlab, y = ylab) +
      scale_y_continuous(limits=c(min(loadings[,dim2],y)-sd(loadings[,dim2]),
                                  max(loadings[,dim2],y)+sd(loadings[,dim2])))+
      scale_x_continuous(limits=c(min(loadings[,dim1],x)-sd(loadings[,dim1]),
                                  max(loadings[,dim1],x)+sd(loadings[,dim1])))
  }
  
  
  else {
    ggplot() +
      geom_point(data = pfs, aes(x = x, y = y, col = "Partial Factor Scores"),size=0.5) +
      geom_text(data = pfs, aes(x = x, y = y, label = rownames(pfs)), hjust = 0, nudge_x = 0.05,size=2) +
      geom_point(data = loadings, aes(x=loadings[,dim1], y=loadings[,dim2], col = "Loadings"),size=0.5) +
      geom_text(data = loadings, aes(x=loadings[,dim1], y=loadings[,dim2], label = rownames(loadings)), hjust = 0, nudge_x = 0.05, size=2) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0)  +
      scale_color_manual(name= "Legend", breaks=c("Partial Factor Scores", "Loadings"), values = c("red", "black")) +
      theme(text = element_text(size=5))+
      labs(title = title, x = xlab, y = ylab)
  }
}

#' @import ggplot2
#' @title Plot Method for object of class \code{"mfa"}
#' @description Plots graphs of factor scores (the compromise of the tables), of partial factor scores, or of variable loadings given two dimensions
#' @param mfa object to be plotted
#' @param dim1 numeric value that determines the 1st PC the user chooses
#' @param dim2 numeric value that determines the 2nd PC the user chooses
#' @param type numeric value that:
#'                if it is 1, then plot the factor scores
#'                if it is 2, then plot the partial factor scores
#'                if it is 3, then plot the loadings
#'                if it is 4, then plot the partial factor scores and variable loadings
#' @param scale_x tunning parameter on the x axis to adjust the scale of the plot
#' @param scale_y tunning parameter on the y axis to adjust the scale of the plot
#' @param text lables for variables
#' @param cat categories to group the observations
#' @param var the vector indicating which subtables user wants to see
#' @return Given the choice of type, return graphs of factor scores (the compromise of the tables), of partial factor scores, or of variable loadings
#' @export
#' @examples
#' cat <- list("New Zealand" = 1:4, "Canada" = 5:8, "France" = 9:12)
#' dftext <- c("Cat Pee", "Passion Fruit", "Green Pepper",
#'             "Mineral", "Smoky", "Cirtrus", "Cat Pee",
#'             "Passion Fruit", "Green Pepper", "Mineral",
#'             "Tropical", "Leafy", "Cat Pee", "Passion Fruit", "Green Pepper",
#'             "Mineral", "Grassy", "Flinty", "Cat Pee", "Passion Fruit", "Green Pepper",
#'             "Mineral", "Leafy", "Cat Pee", "Passion Fruit", "Green Pepper", "Mineral",
#'             "Vegetal", "Hay", "Cat Pee", "Passion Fruit", "Green Pepper", "Mineral",
#'             "Melon", "Cat Pee", "Passion Fruit", "Green Pepper", "Mineral", "Cat Pee",
#'             "Passion Fruit", "Green Pepper", "Mineral", "Grass", "Smoky", "Cat Pee",
#'             "Passion Fruit", "Green Pepper", "Mineral", "Peach", "Cat Pee",
#'             "Passion Fruit", "Green Pepper", "Mineral")
#' plot(mfa = wine, cat = cat, dim1 = 1, dim2 = 2, type = 1, text = dftext)
#' plot(mfa = wine, cat = cat, dim1 = 1, dim2 = 2, type = 2, text = dftext)
#' plot(mfa = wine, dim1 = 1, dim2 = 2, type = 3)
plot.mfa <- function(mfa, dim1 = 1, dim2 = 2, type = 1, scale_x = 1, scale_y = 1, text = NULL, cat = NULL, var = NULL) {
  # check argument
  check_dim(mfa, dim1, dim2)
  check_type(type)
  check_text(mfa, text)
  check_cat(mfa, cat)
  check_var(mfa, var)
  check_scale(scale_x, scale_y)
  if (is.null(var)) {
    var <- 1:length(mfa$sets)
  }
  
  # Plot factor scores
  if (type == "1") {
    plot_fs(mfa, dim1, dim2, cat)
  }
  
  # Plot partial factor scores
  else if (type == "2") {
    if (length(var) == 1) {
      plot_pfs(mfa, dim1, dim2, var, cat)
    }
    else {
      plots <- list()
      for (i in 1:length(var)) {
        plots[[i]] <- plot_pfs(mfa, dim1, dim2, var[i], cat)
      }
      do.call(gridExtra::grid.arrange, c(plots, ncol=sqrt(length(mfa$sets))))
    }
  }
  
  else if (type == "3") {
    if (length(var) == 1) {
      plot_loadings(mfa, dim1, dim2, var, scale_x, scale_y, text)
    }
    else {
      plots <- list()
      for (i in 1:length(var)) {
        plots[[i]] <- plot_loadings(mfa, dim1, dim2, var[i], scale_x, scale_y, text)
      }
      do.call(gridExtra::grid.arrange, c(plots, ncol=sqrt(length(mfa$sets))))
    }
    
  }
  
  else if (type == "4") {
    if (length(var) == 1) {
      plot_pfsl(mfa, dim1, dim2, var, scale_x, scale_y, text, cat)
    }
    else {
      plots <- list()
      for (i in 1:length(var)) {
        plots[[i]] <- plot_pfsl(mfa, dim1, dim2, var[i], scale_x, scale_y, text, cat)
      }
      do.call(gridExtra::grid.arrange, c(plots, ncol=sqrt(length(mfa$sets))))
    }
  }
}
