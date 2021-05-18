########################################################################
######################## BerG-ARMA process #############################
# Code author: Lucas Sales - lucasofsales@gmail.com or lucasofs@ime.usp.br
# Institution: University of São Paulo, SP, Brazil
########################################################################


# The functions qberg and rberg are availabe in a R github package and can be install in 
# two steps:
# first, run the code require("devtools")
# second, run the code devtools::install_github("rdmatheus/bergrm") 
# The bergrm package was developed by Rodrigo Medeiros, Ph.D. student at University of SP


#Generate a BerG-ARMA(p,q) process - Please, be atention to the constrains in the parameters mu and nu.
#The function return a time series object with size n.
#This generation function assumes that the covariate has a senoide form

rbergARMA <- function(n,alpha,beta,phi=NA,theta=NA,
                      nu,freq=12,link="log") {
  qberg <- function(p, mu, nu, lower.tail = TRUE){
    if ((any(p < 0)) || (any(p > 1)))
      stop("p must be in the unit interval: (0, 1)")
    if ((any(mu <= 0)) || (any(nu <= 0)))
      stop("The parameters must be positives")
    if (any(nu < abs(mu - 1)))
      warning("Constraints are not satisfied")
    
    if(lower.tail == FALSE)
      p <- 1 - p
    
    p0 <- (1 - mu + nu)/(1 + mu + nu)
    
    ifelse(length(p) > 1, p.star <- p[p > p0], p.star <- p)
    ifelse(length(mu) > 1, mu.star <- mu[p > p0], mu.star <- mu)
    ifelse(length(nu) > 1, nu.star <- nu[p > p0], nu.star <- nu)
    
    q <- ceiling(
      round(log((1 - p.star) * (1 + mu.star + nu.star) / (2 * mu.star)) /
              log((mu.star + nu.star - 1) / (mu.star + nu.star + 1)), 2)
    )
    
    quanti <- c(rep(0, sum(p <= p0)), q)
    index <- c(which(p <= p0), which(p > p0))
    
    return(quanti[sort(index, index.return = TRUE)$ix])
  }
  
  rberg <- function(n, mu, nu){
    if ((any(mu<=0)) || (any(nu<=0)))
      stop("The parameters must be positives")
    if (any(nu < abs(mu - 1)))
      warning("Constraints are not satisfied")
    
    u <- stats::runif(n)
    return(qberg(u, mu, nu))
  }
  
  ar<-NA
  ma<-NA
  const=0.1
  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }
  
  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("log", "identity")))
  {  
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"log\"  and \"identity\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  ## BerG - ARMA model(p,q)
  if(any(is.na(phi)==F) && any(is.na(theta)==F))
  {
    print("BerG-ARMA model")
    
    p <- max(ar)
    q <- max(ma)
    m <- 100 # this part is add for the burnout in the future
    maxx=max(p,q)
    
    ynew <-rep(alpha,(n+m)) #The initial value of y is the parameter alpha
    mu <- linkinv(ynew)   # The initial value of mu is the g-1(alpha)
    
    error<-rep(0,n+m) # E(error)=0 
    eta<- y <- NULL
    X <- matrix(sin(2*pi*(1:(n+m))/12)) #assuming the covariate as the senoide form
    
    for(i in (maxx+1):(n+m))
    {
      eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta))) + (theta%*%error[i-ma]) 
      mu[i] <- linkinv(eta[i])
      y[i]=rberg(1,mu=mu[i],nu=nu)
      if(y[i]!=0){ynew[i] <- linkfun(y[i])
      } else {ynew[i] <- linkfun(const)}
      error[i]<- ynew[i]-eta[i]   
    }
    if (any(nu < abs(mu - 1)))
      warning("Constraints are not satisfied")
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } # end of ARMA part 
  
  # BerG-AR model
  if(any(is.na(phi)==F) && any(is.na(theta)==T))
  {
    print("BerG-AR model")    
    
    p <- max(ar)
    m <- 100  #this part is add for the burnout in the future
    
    ynew <-rep(alpha,(n+m)) #The initial value of y is the parameter alpha
    mu <- linkinv(ynew)   # The initial value of mu is the g-1(alpha)
    
    eta <- y <- NULL
    X <- cbind(sin(2*pi*(1:(n+m))/12))
    for(i in (p+1):(n+m))
    {
      eta[i]  <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) ))
      mu[i]   <- linkinv(eta[i])
      y[i]=rberg(1,mu=mu[i],nu=nu)
      if(y[i]!=0){ynew[i] <- linkfun(y[i])
      } else {ynew[i] <- linkfun(const)}
    }
    if (any(nu < abs(mu - 1)))
     warning("Constraints are not satisfied")
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } # end of AR part
  
  # MA model
  if(any(is.na(phi)==T) && any(is.na(theta)==F))
  {
    print("MA model")    
    
    q <- max(ma)
    m <- 100 # this part is add for the burnout in the future
    
    ynew <-rep(alpha,(n+m)) #The initial value of y is the parameter alpha
    mu <- linkinv(ynew)   # The initial value of mu is the g-1(alpha)
    
    eta <- y <- error <- rep(0,n+m) # E(error)=0 
    X <- cbind(sin(2*pi*(1:(n+m))/12))
    
    for(i in (q+1):(n+m))
    {
      eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
      mu[i]   <- linkinv(eta[i])
      y[i]=rberg(1,mu=mu[i],nu=nu)
      if(y[i]!=0){ynew[i] <- linkfun(y[i])
      } else {ynew[i] <- linkfun(const)}
      error[i]<- ynew[i]-eta[i]   
      
    }
    if (any(nu < abs(mu - 1)))
     warning("Constraints are not satisfied")
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } # end of MA part  
  
}

class(rbergARMA(100,phi=0.05,beta=0.1,theta=0.05,alpha=0.1,nu=0.5))
