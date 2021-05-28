bergARMA<- function(beta,phi,theta,nu, y, X, link = "log"){
  require(bergrm)
  #Link function: Useful for the conditional Log-Likelihood function, Conditional score and Jacobian matrix
  g <- function(link){
    
    switch(link,
           
           identity = {
             fun <- function(theta) theta
             inv <- function(eta) eta
             deriv. <- function(theta) rep.int(1, length(theta))
             deriv.. <- function(theta) rep.int(0, length(theta))
             valideta <- function(eta) TRUE
           },
           
           log = {
             fun <- function(theta) log(theta)
             inv <- function(eta) pmax(exp(eta), .Machine$double.eps)
             deriv. <- function(theta) 1 / theta
             deriv.. <- function(theta) -1 / (theta ^ 2)
             valideta <- function(eta) TRUE
           },
           
           sqrt = {
             fun <- function(theta) sqrt(theta)
             inv <- function(eta) eta^2
             deriv. <- function(theta) 1 / (2 * sqrt(theta))
             deriv.. <- function(theta) -1 / (4 * (theta ^ (3 / 2)))
             valideta <- function(eta) all(is.finite(eta)) && all(eta > 0)
           },
           
           stop(gettextf("link %s not available", sQuote(link)), domain = NA))
    
    environment(fun) <- environment(inv) <- environment(deriv.) <-
      environment(deriv..) <- environment(valideta) <- asNamespace("stats")
    
    structure(list(fun = fun, inv = inv, deriv. = deriv.,
                   deriv.. = deriv.., valideta = valideta,
                   name = link), class = "link-sdlrm")
  }
  
##### ARMA model 
if(phi != 0 && theta != 0){
  print("ARMA model")
  
  #Necessary quantities for the CLL function, score and jacobian functions
  par<- c(beta,phi,theta,nu)  
  n=length(y)
  p=length(phi);q=length(theta)
  m= max(p,q)
  X <- as.matrix(X)
  const=0.1
  #Conditional LL function
  cond_ll_berg <- function(par, y, X,n,p,q, link = "log"){
    # Creating the link functions
    g1 <- g(link)$fun
    g1.inv <- g(link)$inv
    
    # Necessary quantities
    const=0.1
    X <- as.matrix(X)
    ynew <- g1.inv(y); 
    eta=error=rep(0,length(y))
    m=max(p,q)
    
    beta<- par[1:(ncol(X))]
    phi<-par[(ncol(X)+1):(ncol(X)+p)]
    theta<-par[(ncol(X)+p+1):(ncol(X)+p+q)]
    nu<-par[length(par)]
    #Creating the mu_t; obs: here, I put an threshold when ynew = 0 
    for(i in (m+1):n)
    {
      eta[i] <- X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-(1:p)]-X[i-(1:p),]%*%as.matrix(beta))) + (theta%*%error[i-(1:q)])
      if(y[i]!=0){ynew[i] <- g1(y[i])
      } else {ynew[i] <- g1(const)}
      error[i]<- ynew[i]-eta[i]   
    }
    mu <- g1.inv(eta[(m+1):n])
    ycond <- y[(m+1):n]
    
    
    l0 <- as.numeric(log(1 - mu[ycond == 0] + nu) - log(1 + mu[ycond == 0]+nu))
    l  <- as.numeric(log(4 * mu[ycond > 0]) + (ycond[ycond > 0] - 1)*log(mu[ycond > 0]+nu - 1) -
                       (ycond[ycond > 0] + 1)*log(mu[ycond > 0]+nu + 1))
    #Return the conditional log-likelihood; obs: in the optmization functions I need to put -cond_ll_berg
    return(sum(c(l0,l)))
  }
  
  #Conditional Score function
  cond_U_berg <- function(par,y,X,n,p,q,link="log"){
    # Creating the link functions
    g1 <- g(link)$fun
    g1.inv <- g(link)$inv
    g1. <- g(link)$deriv. #derivative of the link function
    
    # Necessary quantities
    const=0.1
    X <- as.matrix(X)
    ynew <- g1.inv(y); 
    eta=error=rep(0,length(y))
    m=max(p,q)
    
    beta<- par[1:(ncol(X))]
    phi<-par[(ncol(X)+1):(ncol(X)+p)]
    theta<-par[(ncol(X)+p+1):(ncol(X)+p+q)]
    nu<-par[length(par)]
    delta <- (as.numeric(y==0))
    for(i in (m+1):n)
    {
      eta[i] <- X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-(1:p)]-X[i-(1:p),]%*%as.matrix(beta))) + (theta%*%error[i-(1:q)])
      if(y[i]!=0){ynew[i] <- g1(y[i])
      } else {ynew[i] <- g1(const)}
      error[i]<- ynew[i]-eta[i]   
    }
    mu <- g1.inv(eta)
    
    
    
    B_aux <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
    for(i in 1:(n-m))
    {
      for(j in 1:length(beta))
        B_aux[i,j] <- X[i+m,j]-(phi%*%X[i+m-(1:p),j])
    }
    
    P_aux <- matrix(rep(NA,(n-m)*p),ncol=p)
    for(i in 1:(n-m))
    {
      P_aux[i,] <- ynew[i+m-(1:p)] - X[i+m-(1:p),]%*%as.matrix(beta)
    }
    
    T_aux <- matrix(rep(NA,(n-m)*q),ncol=q)
    for(i in 1:(n-m))
    {
      T_aux[i,] <- error[i+m-(1:q)]
    }
    
    V_aux <- -2*delta*(1+nu)/((1-mu+nu)*(1+mu+nu)) + (1-delta)*(1/mu + 2*(y-mu-nu)/(mu+nu-1)*(mu+nu+1))
    
    
    ###Initializing the derivates in zero
    #deta.dalpha <- matrix(0, ncol=1,nrow=n)
    deta.dbeta <- matrix(0, ncol=length(beta),nrow=n)
    deta.dphi <- matrix(0, ncol=p,nrow=n)
    deta.dtheta <- matrix(0, ncol=q,nrow=n)
    
    for(i in (m+1):n)
    {
      #deta.dalpha[i,] <- 1 -  theta%*%deta.dalpha[i-(1:q),]
      deta.dbeta[i,]<- B_aux[(i-m),] - theta%*%deta.dbeta[i-(1:q),]
      deta.dphi[i,] <- P_aux[(i-m),] - theta%*%deta.dphi[i-(1:q),] 
      deta.dtheta[i,]<- T_aux[(i-m),] - theta%*%deta.dtheta[i-(1:q),]
    }
    ### Derivate with relation to nu
    deta.dnu <- 2*delta*mu/((1-mu+nu)*(1+mu+nu)) + (1-delta)*( 2*(y-mu-nu)/(mu+nu-1)*(mu+nu+1))
    #Matrices use in the score vector (matrix Q = matrix T of the paper)
    #A <- deta.dalpha[(m+1):n]
    B <- deta.dbeta[(m+1):n,]
    P <- deta.dphi[(m+1):n,]
    Q <- deta.dtheta[(m+1):n,]
    D <- diag(as.numeric(1/g1.(mu[(m+1):n])))
    V <- V_aux[(m+1):n]
    
    #Ua <- t(A)%*%D%*%V
    Ub <- t(B)%*%D%*%V
    Up <- t(P)%*%D%*%V
    Ut <- t(Q)%*%D%*%V
    Un <- sum(deta.dnu[(m+1): n])
    
    U<- c(Ub,Up,Ut,Un)   
    return(U)
  }
  
  
  ll <- function(par) -cond_ll_berg(par, y, X,n,p,q, link = "log")
  U <- function(par) -cond_U_berg(par,y,X,n,p,q,link="log")
  
  opt <- nlminb(par, ll, gradient = NULL, hessian = NULL, scale = 1, control = list(),
                lower = c(rep(-Inf,ncol(X)),rep(-1,p+q),0.001), 
                upper = c(rep(Inf,ncol(X)),rep(1,p+q),Inf))
  
  
  aux <- c()
  aux$serie <- y
  aux$par <- opt$par
  names(aux$par)<-c(paste("beta",0:(ncol(X)-1)),paste("phi",1:p),paste("theta",1:q),"nu")
  aux$loglik <- opt$objective  
  aux$AIC <- 2*aux$loglik+2*(ncol(X)+p+q+1)
  aux$BIC <- 2*aux$loglik+log(n)*(ncol(X)+p+q+1)
  
  ## residuals
  
  #Will be used in the residuals 
  mu_residuals <- function(par,y,X){  
    # Creating the link functions
    g1 <- g(link)$fun
    g1.inv <- g(link)$inv
    beta.hat<- par[1:(ncol(X))]
    phi.hat<-par[(ncol(X)+1):(ncol(X)+p)]
    theta.hat<-par[(ncol(X)+p+1):(ncol(X)+p+q)]
    nu.hat<-par[length(aux$par)]
    X <- as.matrix(X)
    ynew <- g1.inv(y); 
    eta<-error<-rep(0,length(y))
    for(i in (m+1):n)
    {
      eta[i] <- X[i,]%*%as.matrix(beta.hat) + (phi.hat%*%(ynew[i-(1:p)]-X[i-(1:p),]%*%as.matrix(beta.hat))) + (theta.hat%*%error[i-(1:q)])
      if(y[i]!=0){ynew[i] <- g1(y[i])
      } else {ynew[i] <- g1(const)}
      error[i]<- ynew[i]-eta[i]   
    }
    mu.hat <- g1.inv(eta)
    return(mu.hat)
  }
  
  
  ut <-runif(n-m)
  mu.hat=mu_residuals(aux$par,y,X)    
  aux$mu <- mu.hat
  p_berg <- d_berg <- NULL
  
  for (i in (m+1):n) {
    d_berg[i] <-dberg(y[i],mu.hat[i],aux$par[length(aux$par)])
    p_berg[i] <- pberg(y[i]-1,mu.hat[i],aux$par[length(aux$par)])
  }
  aux$resid <- qnorm( p_berg[(m+1):n] + ut*d_berg[(m+1):n])
  
  ## covariance matrix
  berg_cov <- function(par,y,X,n,p,q,link="log"){
    # Creating the link functions
    g1 <- g(link)$fun
    g1.inv <- g(link)$inv
    g1. <- g(link)$deriv. #derivative of the link function
    
    # Necessary quantities
    const=0.1
    X <- as.matrix(X)
    ynew <- g1.inv(y); 
    eta=error=rep(0,length(y))
    m=max(p,q)
    
    beta<- par[1:(ncol(X))]
    phi<-par[(ncol(X)+1):(ncol(X)+p)]
    theta<-par[(ncol(X)+p+1):(ncol(X)+p+q)]
    nu<-par[length(par)]
    delta <- (as.numeric(y==0))
    for(i in (m+1):n)
    {
      eta[i] <- X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-(1:p)]-X[i-(1:p),]%*%as.matrix(beta))) + (theta%*%error[i-(1:q)])
      if(y[i]!=0){ynew[i] <- g1(y[i])
      } else {ynew[i] <- g1(const)}
      error[i]<- ynew[i]-eta[i]   
    }
    mu <- g1.inv(eta)
    
    
    
    B_aux <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
    for(i in 1:(n-m))
    {
      for(j in 1:length(beta))
        B_aux[i,j] <- X[i+m,j]-(phi%*%X[i+m-(1:p),j])
    }
    
    P_aux <- matrix(rep(NA,(n-m)*p),ncol=p)
    for(i in 1:(n-m))
    {
      P_aux[i,] <- ynew[i+m-(1:p)] - X[i+m-(1:p),]%*%as.matrix(beta)
    }
    
    T_aux <- matrix(rep(NA,(n-m)*q),ncol=q)
    for(i in 1:(n-m))
    {
      T_aux[i,] <- error[i+m-(1:q)]
    }
    
    
    w <- 4*mu*(1+nu)/((1-mu+nu)*(1+mu+nu)^3) + 2/(mu*(1+mu+nu)) + 4*mu/((mu+nu-1)*(mu+nu+1)^3)
    #Kvv is the sum(v1 in m+1 to n) and v2 is useful for the cov between nu and the parameters
    v1 <- 4*mu/((mu+nu-1)*(mu+nu+1)^3) + 4*mu*(nu+1)/((1-mu+nu)*(1+mu+nu)^3)
    v2 <- (2*(nu+1)^2 + 2*mu^2)/((1-mu+nu)*(1+mu+nu)^3) - 4*mu*(nu+1)/((mu+nu-1)*(mu+nu+1)^3) 
    ###Initializing the derivates in zero
    #deta.dalpha <- matrix(0, ncol=1,nrow=n)
    deta.dbeta <- matrix(0, ncol=length(beta),nrow=n)
    deta.dphi <- matrix(0, ncol=p,nrow=n)
    deta.dtheta <- matrix(0, ncol=q,nrow=n)
    dmu.deta <- as.numeric(1/g1.(mu)) 
    for(i in (m+1):n)
    {
      #deta.dalpha[i,] <- 1 -  theta%*%deta.dalpha[i-(1:q),]
      deta.dbeta[i,]<- B_aux[(i-m),] - theta%*%deta.dbeta[i-(1:q),]
      deta.dphi[i,] <- P_aux[(i-m),] - theta%*%deta.dphi[i-(1:q),] 
      deta.dtheta[i,]<- T_aux[(i-m),] - theta%*%deta.dtheta[i-(1:q),]
    }
    ### Derivate with relation to nu
    deta.dnu <- 2*delta*mu/((1-mu+nu)*(1+mu+nu)) + (1-delta)*( 2*(y-mu-nu)/(mu+nu-1)*(mu+nu+1))
    #Matrices use in the score vector (matrix Q = matrix T of the paper)
    #A <- deta.dalpha[(m+1):n]
    B <- deta.dbeta[(m+1):n,]
    P <- deta.dphi[(m+1):n,]
    Q <- deta.dtheta[(m+1):n,]
    D <- diag(as.numeric(1/g1.(mu[(m+1):n])))
    W <- diag(w[(m+1):n])
    V <- matrix(v2[(m+1):n],ncol=(n-m))
    
    Kbb <- t(B)%*%W%*%D%*%D%*%B
    Kbp <- t(B)%*%W%*%D%*%D%*%P
    Kbt <- t(B)%*%W%*%D%*%D%*%Q
    Kbv <- t(V%*%D%*%B)
    Kpp <- t(P)%*%W%*%D%*%D%*%P
    Kpb <- t(Kbp)
    Kpt <- t(P)%*%W%*%D%*%D%*%Q
    Kpv <- t(V%*%D%*%P)
    Ktt <- t(Q)%*%W%*%D%*%D%*%Q 
    Ktb <- t(Kbt)
    Ktp <- t(Kpt)
    Ktv <- t(V%*%D%*%Q)
    Kvv <- sum(v1[(m+1):n])
    Kvb <- t(Kbv)
    Kvp <- t(Kpv)
    Kvt <- t(Ktv)
    
  K<- rbind(  cbind(Kbb,Kbp,Kbt,Kbv),cbind(Kpb,Kpp,Kpt,Kpv)
    ,cbind(Ktb,Ktp,Ktt,Ktv),cbind(Kvb,Kvp,Kvt,Kvv))
    
  return(solve(K))
  }
  
  aux$cov<-berg_cov(aux$par,y,X,n,p,q,link="log")
  
} #Final of ARMA Model
  
  ##### AR model 
  if(phi != 0 && theta == 0){
    print("AR model")
    
    #Necessary quantities for the CLL function, score and jacobian functions
    par<- c(beta,phi,nu)  
    n=length(y)
    p=length(phi);q=0
    m= max(p,q)
    X <- as.matrix(X)
    const=0.1
    #Conditional LL function
    cond_ll_berg <- function(par, y, X,n,p,q, link = "log"){
      # Creating the link functions
      g1 <- g(link)$fun
      g1.inv <- g(link)$inv
      
      # Necessary quantities
      const=0.1
      X <- as.matrix(X)
      ynew <- g1.inv(y); 
      eta=error=rep(0,length(y))
      m=max(p,q)
      
      beta<- par[1:(ncol(X))]
      phi<-par[(ncol(X)+1):(ncol(X)+p)]
      nu<-par[length(par)]
      #Creating the mu_t; obs: here, I put an threshold when ynew = 0 
      for(i in (m+1):n)
      {
        eta[i] <- X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-(1:p)]-X[i-(1:p),]%*%as.matrix(beta))) 
        if(y[i]!=0){ynew[i] <- g1(y[i])
        } else {ynew[i] <- g1(const)}
        error[i]<- ynew[i]-eta[i]   
      }
      mu <- g1.inv(eta[(m+1):n])
      ycond <- y[(m+1):n]
      
      
      l0 <- as.numeric(log(1 - mu[ycond == 0] + nu) - log(1 + mu[ycond == 0]+nu))
      l  <- as.numeric(log(4 * mu[ycond > 0]) + (ycond[ycond > 0] - 1)*log(mu[ycond > 0]+nu - 1) -
                         (ycond[ycond > 0] + 1)*log(mu[ycond > 0]+nu + 1))
      #Return the conditional log-likelihood; obs: in the optmization functions I need to put -cond_ll_berg
      return(sum(c(l0,l)))
    }
    
    #Conditional Score function
    cond_U_berg <- function(par,y,X,n,p,q,link="log"){
      # Creating the link functions
      g1 <- g(link)$fun
      g1.inv <- g(link)$inv
      g1. <- g(link)$deriv. #derivative of the link function
      
      # Necessary quantities
      const=0.1
      X <- as.matrix(X)
      ynew <- g1.inv(y); 
      eta=error=rep(0,length(y))
      m=max(p,q)
      
      beta<- par[1:(ncol(X))]
      phi<-par[(ncol(X)+1):(ncol(X)+p)]
      nu<-par[length(par)]
      delta <- (as.numeric(y==0))
      
      for(i in (m+1):n)
      {
        eta[i] <- X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-(1:p)]-X[i-(1:p),]%*%as.matrix(beta))) 
        if(y[i]!=0){ynew[i] <- g1(y[i])
        } else {ynew[i] <- g1(const)}
        error[i]<- ynew[i]-eta[i]   
      }
      mu <- g1.inv(eta)
      
      
      
      B_aux <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
      for(i in 1:(n-m))
      {
        for(j in 1:length(beta))
          B_aux[i,j] <- X[i+m,j]-(phi%*%X[i+m-(1:p),j])
      }
      
      P_aux <- matrix(rep(NA,(n-m)*p),ncol=p)
      for(i in 1:(n-m))
      {
        P_aux[i,] <- ynew[i+m-(1:p)] - X[i+m-(1:p),]%*%as.matrix(beta)
      }
      
      
      
      V_aux <- -2*delta*(1+nu)/((1-mu+nu)*(1+mu+nu)) + (1-delta)*(1/mu + 2*(y-mu-nu)/(mu+nu-1)*(mu+nu+1))
      
      
      ###Initializing the derivates in zero
      #deta.dalpha <- matrix(0, ncol=1,nrow=n)
      deta.dbeta <- matrix(0, ncol=length(beta),nrow=n)
      deta.dphi <- matrix(0, ncol=p,nrow=n)
      
      for(i in (m+1):n)
      {
        #deta.dalpha[i,] <- 1 
        deta.dbeta[i,]<- B_aux[(i-m),]
        deta.dphi[i,] <- P_aux[(i-m),] 
      }
      ### Derivate with relation to nu
      deta.dnu <- 2*delta*mu/((1-mu+nu)*(1+mu+nu)) + (1-delta)*( 2*(y-mu-nu)/(mu+nu-1)*(mu+nu+1))
      #Matrices use in the score vector (matrix Q = matrix T of the paper)
      #A <- deta.dalpha[(m+1):n]
      B <- deta.dbeta[(m+1):n,]
      P <- deta.dphi[(m+1):n,]
      D <- diag(as.numeric(1/g1.(mu[(m+1):n])))
      V <- V_aux[(m+1):n]
      
      #Ua <- t(A)%*%D%*%V
      Ub <- t(B)%*%D%*%V
      Up <- t(P)%*%D%*%V
      Un <- sum(deta.dnu[(m+1): n])
      
      U<- c(Ub,Up,Un)   
      return(U)
    }
    
    
    ll <- function(par) -cond_ll_berg(par, y, X,n,p,q, link = "log")
    U <- function(par) -cond_U_berg(par,y,X,n,p,q,link="log")
    
    opt <- nlminb(par, ll, gradient = NULL, hessian = NULL, scale = 1, control = list(),
                  lower = c(rep(-Inf,ncol(X)),rep(-1,p+q),0.001), 
                  upper = c(rep(Inf,ncol(X)),rep(1,p+q),Inf))
    
    
    aux <- c()
    aux$serie <- y
    aux$par <- opt$par
    names(aux$par)<-c(paste("beta",0:(ncol(X)-1)),paste("phi",1:p),"nu")
    aux$loglik <- opt$objective  
    aux$AIC <- 2*aux$loglik+2*(ncol(X)+p+q+1)
    aux$BIC <- 2*aux$loglik+log(n)*(ncol(X)+p+q+1)
    
    ## residuals
    
    #Will be used in the residuals 
    mu_residuals <- function(par,y,X){  
      # Creating the link functions
      g1 <- g(link)$fun
      g1.inv <- g(link)$inv
      beta.hat<- par[1:(ncol(X))]
      phi.hat<-par[(ncol(X)+1):(ncol(X)+p)]
      nu.hat<-par[length(aux$par)]
      X <- as.matrix(X)
      ynew <- g1.inv(y); 
      eta<-error<-rep(0,length(y))
      for(i in (m+1):n)
      {
        eta[i] <- X[i,]%*%as.matrix(beta.hat) + (phi.hat%*%(ynew[i-(1:p)]-X[i-(1:p),]%*%as.matrix(beta.hat)))
        if(y[i]!=0){ynew[i] <- g1(y[i])
        } else {ynew[i] <- g1(const)}
        error[i]<- ynew[i]-eta[i]   
      }
      mu.hat <- g1.inv(eta)
      return(mu.hat)
    }
    
    ut <-runif(n-m)
    mu.hat=mu_residuals(aux$par,y,X)    
    aux$mu <- mu.hat
    p_berg <- d_berg <- NULL
    for (i in (m+1):n) {
      d_berg[i] <-dberg(y[i],mu.hat[i],aux$par[length(aux$par)])
      p_berg[i] <- pberg(y[i]-1,mu.hat[i],aux$par[length(aux$par)])
    }
    
    aux$resid <- qnorm( p_berg[(m+1):n] + ut*d_berg[(m+1):n])
    
      ## covariance matrix
  berg_cov <- function(par,y,X,n,p,q,link="log"){
 # Creating the link functions
      g1 <- g(link)$fun
      g1.inv <- g(link)$inv
      g1. <- g(link)$deriv. #derivative of the link function
      
      # Necessary quantities
      const=0.1
      X <- as.matrix(X)
      ynew <- g1.inv(y); 
      eta=error=rep(0,length(y))
      m=max(p,q)
      
      beta<- par[1:(ncol(X))]
      phi<-par[(ncol(X)+1):(ncol(X)+p)]
      nu<-par[length(par)]
      delta <- (as.numeric(y==0))
      
      for(i in (m+1):n)
      {
        eta[i] <- X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-(1:p)]-X[i-(1:p),]%*%as.matrix(beta))) 
        if(y[i]!=0){ynew[i] <- g1(y[i])
        } else {ynew[i] <- g1(const)}
        error[i]<- ynew[i]-eta[i]   
      }
      mu <- g1.inv(eta)
      
      
      
      B_aux <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
      for(i in 1:(n-m))
      {
        for(j in 1:length(beta))
          B_aux[i,j] <- X[i+m,j]-(phi%*%X[i+m-(1:p),j])
      }
      
      P_aux <- matrix(rep(NA,(n-m)*p),ncol=p)
      for(i in 1:(n-m))
      {
        P_aux[i,] <- ynew[i+m-(1:p)] - X[i+m-(1:p),]%*%as.matrix(beta)
      }
      
      
      
      ###Initializing the derivates in zero
      #deta.dalpha <- matrix(0, ncol=1,nrow=n)
      deta.dbeta <- matrix(0, ncol=length(beta),nrow=n)
      deta.dphi <- matrix(0, ncol=p,nrow=n)
      
      for(i in (m+1):n)
      {
        #deta.dalpha[i,] <- 1 
        deta.dbeta[i,]<- B_aux[(i-m),]
        deta.dphi[i,] <- P_aux[(i-m),] 
      }
      ### Derivate with relation to nu
      deta.dnu <- 2*delta*mu/((1-mu+nu)*(1+mu+nu)) + (1-delta)*( 2*(y-mu-nu)/(mu+nu-1)*(mu+nu+1))
      #Matrices use in the score vector (matrix Q = matrix T of the paper)
      #A <- deta.dalpha[(m+1):n]
      B <- deta.dbeta[(m+1):n,]
      P <- deta.dphi[(m+1):n,]
      D <- diag(as.numeric(1/g1.(mu[(m+1):n])))
      
    
    w <- 4*mu*(1+nu)/((1-mu+nu)*(1+mu+nu)^3) + 2/(mu*(1+mu+nu)) + 4*mu/((mu+nu-1)*(mu+nu+1)^3)
    #Kvv is the sum(v1 in m+1 to n) and v2 is useful for the cov between nu and the parameters
    v1 <- 4*mu/((mu+nu-1)*(mu+nu+1)^3) + 4*mu*(nu+1)/((1-mu+nu)*(1+mu+nu)^3)
    v2 <- (2*(nu+1)^2 + 2*mu^2)/((1-mu+nu)*(1+mu+nu)^3) - 4*mu*(nu+1)/((mu+nu-1)*(mu+nu+1)^3) 
    V <- matrix(v2[(m+1):n],ncol=(n-m))
    W <- diag(w[(m+1):n])
    
    Kbb <- t(B)%*%W%*%D%*%D%*%B
    Kbp <- t(B)%*%W%*%D%*%D%*%P
    Kbv <- t(V%*%D%*%B)
    Kpp <- t(P)%*%W%*%D%*%D%*%P
    Kpb <- t(Kbp)
    Kpv <- t(V%*%D%*%P)
    Kvv <- sum(v1[(m+1):n])
    Kvb <- t(Kbv)
    Kvp <- t(Kpv)
   
    
  K<- rbind(  cbind(Kbb,Kbp,Kbv),cbind(Kpb,Kpp,Kpv),
    cbind(Kvb,Kvp,Kvv))
    
  return(solve(K))
  }
  
  aux$cov<-berg_cov(aux$par,y,X,n,p,q,link="log")
    
    
  } #Final of AR Model
  
  ##### MA model 
  if(phi == 0 && theta != 0){
    print("MA model")
    
    #Necessary quantities for the CLL function, score and jacobian functions
    par<- c(beta,theta,nu)  
    n=length(y)
    p=0;q=length(theta)
    m= max(p,q)
    X <- as.matrix(X)
    const=0.1
    #Conditional LL function
    cond_ll_berg <- function(par, y, X,n,p,q, link = "log"){
      # Creating the link functions
      g1 <- g(link)$fun
      g1.inv <- g(link)$inv
      
      # Necessary quantities
      const=0.1
      X <- as.matrix(X)
      ynew <- g1.inv(y); 
      eta=error=rep(0,length(y))
      m=max(p,q)
      
      
      beta<- par[1:(ncol(X))]
      theta<-par[(ncol(X)+1):(ncol(X)+q)]
      nu<-par[length(par)]
      #Creating the mu_t; obs: here, I put an threshold when ynew = 0 
      for(i in (m+1):n)
      {
        eta[i] <- X[i,]%*%as.matrix(beta)  + (theta%*%error[i-(1:q)])
        if(y[i]!=0){ynew[i] <- g1(y[i])
        } else {ynew[i] <- g1(const)}
        error[i]<- ynew[i]-eta[i]   
      }
      mu <- g1.inv(eta[(m+1):n])
      ycond <- y[(m+1):n]
      
      
      l0 <- as.numeric(log(1 - mu[ycond == 0] + nu) - log(1 + mu[ycond == 0]+nu))
      l  <- as.numeric(log(4 * mu[ycond > 0]) + (ycond[ycond > 0] - 1)*log(mu[ycond > 0]+nu - 1) -
                         (ycond[ycond > 0] + 1)*log(mu[ycond > 0]+nu + 1))
      #Return the conditional log-likelihood; obs: in the optmization functions I need to put -cond_ll_berg
      return(sum(c(l0,l)))
    }
    
    #Conditional Score function
    cond_U_berg <- function(par,y,X,n,p,q,link="log"){
      # Creating the link functions
      g1 <- g(link)$fun
      g1.inv <- g(link)$inv
      g1. <- g(link)$deriv. #derivative of the link function
      
      # Necessary quantities
      const=0.1
      X <- as.matrix(X)
      ynew <- g1.inv(y); 
      eta=error=rep(0,length(y))
      m=max(p,q)
      
      beta<- par[1:(ncol(X))]
      theta<-par[(ncol(X)+1):(ncol(X)+q)]
      nu<-par[length(par)]
      delta <- (as.numeric(y==0))
      
      for(i in (m+1):n)
      {
        eta[i] <- X[i,]%*%as.matrix(beta) + (theta%*%error[i-(1:q)])
        if(y[i]!=0){ynew[i] <- g1(y[i])
        } else {ynew[i] <- g1(const)}
        error[i]<- ynew[i]-eta[i]   
      }
      mu <- g1.inv(eta)
      
      
      
      
      B_aux <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
      for(i in 1:(n-m))
      {
        for(j in 1:length(beta))
          B_aux[i,j] <- X[i+m,j]
      }
      
      
      T_aux <- matrix(rep(NA,(n-m)*q),ncol=q)
      for(i in 1:(n-m))
      {
        T_aux[i,] <- error[i+m-(1:q)]
      }
      
      V_aux <- -2*delta*(1+nu)/((1-mu+nu)*(1+mu+nu)) + (1-delta)*(1/mu + 2*(y-mu-nu)/(mu+nu-1)*(mu+nu+1))
      
      
      ###Initializing the derivates in zero
      #deta.dalpha <- matrix(0, ncol=1,nrow=n)
      deta.dbeta <- matrix(0, ncol=length(beta),nrow=n)
      deta.dphi <- matrix(0, ncol=p,nrow=n)
      deta.dtheta <- matrix(0, ncol=q,nrow=n)
      
      for(i in (m+1):n)
      {
        #deta.dalpha[i,] <- 1 -  theta%*%deta.dalpha[i-(1:q),]
        deta.dbeta[i,]<- B_aux[(i-m),] - theta%*%deta.dbeta[i-(1:q),]
        deta.dtheta[i,]<- T_aux[(i-m),] - theta%*%deta.dtheta[i-(1:q),]
      }
      ### Derivate with relation to nu
      deta.dnu <- 2*delta*mu/((1-mu+nu)*(1+mu+nu)) + (1-delta)*( 2*(y-mu-nu)/(mu+nu-1)*(mu+nu+1))
      #Matrices use in the score vector (matrix Q = matrix T of the paper)
      #A <- deta.dalpha[(m+1):n]
      B <- deta.dbeta[(m+1):n,]
      Q <- deta.dtheta[(m+1):n,]
      D <- diag(as.numeric(1/g1.(mu[(m+1):n])))
      V <- V_aux[(m+1):n]
      
      #Ua <- t(A)%*%D%*%V
      Ub <- t(B)%*%D%*%V
      Ut <- t(Q)%*%D%*%V
      Un <- sum(deta.dnu[(m+1): n])
      
      U<- c(Ub,Ut,Un)   
      return(U)
    }
    
    
    ll <- function(par) -cond_ll_berg(par, y, X,n,p,q, link = "log")
    U <- function(par) -cond_U_berg(par,y,X,n,p,q,link="log")
    
    opt <- nlminb(par, ll, gradient = NULL, hessian = NULL, scale = 1, control = list(),
                  lower = c(rep(-Inf,ncol(X)),rep(-1,p+q),0.001), 
                  upper = c(rep(Inf,ncol(X)),rep(1,p+q),Inf))
    
    
    aux <- c()
    aux$serie <- y
    aux$par <- opt$par
    names(aux$par)<-c(paste("beta",0:(ncol(X)-1)),paste("theta",1:q),"nu")
    aux$loglik <- opt$objective  
    aux$AIC <- 2*aux$loglik+2*(ncol(X)+p+q+1)
    aux$BIC <- 2*aux$loglik+log(n)*(ncol(X)+p+q+1)
    
    ## residuals
    
    #Will be used in the residuals 
    mu_residuals <- function(par,y,X){  
      # Creating the link functions
      g1 <- g(link)$fun
      g1.inv <- g(link)$inv
      beta.hat<- par[1:(ncol(X))]
      theta.hat<-par[(ncol(X)+1):(ncol(X)+q)]
      nu.hat<-par[length(aux$par)]
      X <- as.matrix(X)
      ynew <- g1.inv(y); 
      eta<-error<-rep(0,length(y))
      for(i in (m+1):n)
      {
        eta[i] <- X[i,]%*%as.matrix(beta.hat)  + (theta.hat%*%error[i-(1:q)])
        if(y[i]!=0){ynew[i] <- g1(y[i])
        } else {ynew[i] <- g1(const)}
        error[i]<- ynew[i]-eta[i]   
      }
      mu.hat <- g1.inv(eta)
      return(mu.hat)
    }
    
    ut <-runif(n-m)
    mu.hat=mu_residuals(aux$par,y,X)    
    aux$mu <- mu.hat
    p_berg <- d_berg <- NULL
    for (i in (m+1):n) {
      d_berg[i] <-dberg(y[i],mu.hat[i],aux$par[length(aux$par)])
      p_berg[i] <- pberg(y[i]-1,mu.hat[i],aux$par[length(aux$par)])
    }
    aux$resid <- qnorm( p_berg[(m+1):n] + ut*d_berg[(m+1):n])
    
     berg_cov <- function(par,y,X,n,p,q,link="log"){
# Creating the link functions
      g1 <- g(link)$fun
      g1.inv <- g(link)$inv
      g1. <- g(link)$deriv. #derivative of the link function
      
      # Necessary quantities
      const=0.1
      X <- as.matrix(X)
      ynew <- g1.inv(y); 
      eta=error=rep(0,length(y))
      m=max(p,q)
      
      beta<- par[1:(ncol(X))]
      theta<-par[(ncol(X)+1):(ncol(X)+q)]
      nu<-par[length(par)]
      delta <- (as.numeric(y==0))
      
      for(i in (m+1):n)
      {
        eta[i] <- X[i,]%*%as.matrix(beta) + (theta%*%error[i-(1:q)])
        if(y[i]!=0){ynew[i] <- g1(y[i])
        } else {ynew[i] <- g1(const)}
        error[i]<- ynew[i]-eta[i]   
      }
      mu <- g1.inv(eta)
      
      
      
      
      B_aux <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
      for(i in 1:(n-m))
      {
        for(j in 1:length(beta))
          B_aux[i,j] <- X[i+m,j]
      }
      
      
      T_aux <- matrix(rep(NA,(n-m)*q),ncol=q)
      for(i in 1:(n-m))
      {
        T_aux[i,] <- error[i+m-(1:q)]
      }
      

      
      ###Initializing the derivates in zero
      #deta.dalpha <- matrix(0, ncol=1,nrow=n)
      deta.dbeta <- matrix(0, ncol=length(beta),nrow=n)
      deta.dphi <- matrix(0, ncol=p,nrow=n)
      deta.dtheta <- matrix(0, ncol=q,nrow=n)
      
      for(i in (m+1):n)
      {
        #deta.dalpha[i,] <- 1 -  theta%*%deta.dalpha[i-(1:q),]
        deta.dbeta[i,]<- B_aux[(i-m),] - theta%*%deta.dbeta[i-(1:q),]
        deta.dtheta[i,]<- T_aux[(i-m),] - theta%*%deta.dtheta[i-(1:q),]
      }
      ### Derivate with relation to nu
      deta.dnu <- 2*delta*mu/((1-mu+nu)*(1+mu+nu)) + (1-delta)*( 2*(y-mu-nu)/(mu+nu-1)*(mu+nu+1))
      #Matrices use in the score vector (matrix Q = matrix T of the paper)
      #A <- deta.dalpha[(m+1):n]
      B <- deta.dbeta[(m+1):n,]
      Q <- deta.dtheta[(m+1):n,]
      D <- diag(as.numeric(1/g1.(mu[(m+1):n])))
     w <- 4*mu*(1+nu)/((1-mu+nu)*(1+mu+nu)^3) + 2/(mu*(1+mu+nu)) + 4*mu/((mu+nu-1)*(mu+nu+1)^3)
    #Kvv is the sum(v1 in m+1 to n) and v2 is useful for the cov between nu and the parameters
    v1 <- 4*mu/((mu+nu-1)*(mu+nu+1)^3) + 4*mu*(nu+1)/((1-mu+nu)*(1+mu+nu)^3)
    v2 <- (2*(nu+1)^2 + 2*mu^2)/((1-mu+nu)*(1+mu+nu)^3) - 4*mu*(nu+1)/((mu+nu-1)*(mu+nu+1)^3)
    W <- diag(w[(m+1):n])
    V <- matrix(v2[(m+1):n],ncol=(n-m))
    
    Kbb <- t(B)%*%W%*%D%*%D%*%B
    Kbt <- t(B)%*%W%*%D%*%D%*%Q
    Kbv <- t(V%*%D%*%B)
    Ktt <- t(Q)%*%W%*%D%*%D%*%Q 
    Ktb <- t(Kbt)
    Ktv <- t(V%*%D%*%Q)
    Kvv <- sum(v1[(m+1):n])
    Kvb <- t(Kbv)
    Kvt <- t(Ktv)
    
  K<- rbind(  cbind(Kbb,Kbt,Kbv),cbind(Ktb,Ktt,Ktv),cbind(Kvb,Kvt,Kvv))
    
  return(solve(K))
  }
  
  aux$cov<-berg_cov(aux$par,y,X,n,p,q,link="log")
  } #Final of MA Model
  
  return(aux)  
  
  
}


