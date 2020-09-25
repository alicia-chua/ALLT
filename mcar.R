rm(list = ls())

#setwd("C:\\Users\\Alicia\\Desktop")

library(nlme)
library(dplyr)
library(tidyverse)
library(data.table)
library(MASS)


simulation <- function(nsamp, truebeta, minnumvis, sampsize, seed, pctmiss){
  
  
  #nacc <- nacc %>% group_by(naccid) %>% filter(n()==minnumvis) %>% ungroup()
  #nacc.id<-nacc[!duplicated(nacc[c("naccid")]),"naccid"]
  
  # Initialize vectors
  b.lme.rirs.ar1    <-rep(0,nsamp)
  se.rirs.ar1       <-rep(0,nsamp)
  ci.low.rirs.ar1   <-rep(0,nsamp)
  ci.upp.rirs.ar1   <-rep(0,nsamp)
  coverage.rirs.ar1 <-rep(0,nsamp)
  tstat.rirs.ar1    <-rep(0,nsamp)
  
  b.lme.ri.ar1      <-rep(0,nsamp)
  se.ri.ar1         <-rep(0,nsamp)
  ci.low.ri.ar1     <-rep(0,nsamp)
  ci.upp.ri.ar1     <-rep(0,nsamp)
  coverage.ri.ar1   <-rep(0,nsamp)
  tstat.ri.ar1      <-rep(0,nsamp)
  
  b.lme.rirs.un     <-rep(0,nsamp)
  se.rirs.un        <-rep(0,nsamp)
  ci.low.rirs.un    <-rep(0,nsamp)
  ci.upp.rirs.un    <-rep(0,nsamp)
  coverage.rirs.un  <-rep(0,nsamp)
  tstat.rirs.un     <-rep(0,nsamp)
  
  b.lme.ri.un     <-rep(0,nsamp)
  se.ri.un        <-rep(0,nsamp)
  ci.low.ri.un    <-rep(0,nsamp)
  ci.upp.ri.un    <-rep(0,nsamp)
  coverage.ri.un  <-rep(0,nsamp)
  tstat.ri.un     <-rep(0,nsamp)
  
  beta.llt          <-rep(0,nsamp)
  sd.beta.llt       <-rep(0,nsamp)
  se.beta.llt       <-rep(0,nsamp)
  ci.low.llt        <-rep(0,nsamp)
  ci.upp.llt        <-rep(0,nsamp)
  coverage.llt      <-rep(0,nsamp)
  tstat.llt         <-rep(0,nsamp)
  
  set.seed(seed)
  
  skip_to_next <- FALSE 
  
  for (s in 1:nsamp){
    tryCatch(expr = {   
    n <- 4000
    ### average intercept and slope
    beta0 <- 1.0
    beta1 <- truebeta
    ### true autocorrelation
    ar.val <- .8
    ### true error SD, intercept SD, slope SD, and intercept-slope cor
    sigma <- 2
    tau0  <- 2.5
    tau1  <- 2.0
    tau01 <- 0.3
    ### maximum number of possible observations
    m <- minnumvis
    ### simulate number of observations for each individual
    p <- round(runif(n,3,m))
    ### simulate observation moments (assume everybody has 1st obs)
    obs <- unlist(sapply(p, function(x) c(1, sort(sample(2:m, x-1, replace=FALSE)))))
    ### set up data frame
    dat <- data.frame(naccid=rep(1:n, times=p))
    naccid <- unique(dat$naccid)
    group<-rbinom(length(naccid),1,.5)
    id.grp<-data.frame(cbind(naccid,group))
    baseage<-data.frame(scale(rnorm(length(naccid), 72.06, 9.61)))
    baseage<-cbind(baseage,id.grp)
    colnames(baseage)[1] <- "baseage"
    colnames(baseage)[3] <- "group"
    
    
    dat<-merge(baseage, dat, by = "naccid")
    
    tsl<-rnorm(nrow(dat),1,0.2)
    dat<-data.frame(cbind(dat, tsl))
    dat <- dat  %>% group_by(naccid)%>% mutate(vnumber = row_number(),tsl_adj = ifelse(row_number()==1,0,tsl),
                                               timesincefirst = cumsum(tsl_adj)) %>% ungroup()    
    dat <- data.frame(dat)
    colnames(dat)[3] <- "group"
    
    ### simulate (correlated) random effects for intercepts and slopes
    mu  <- c(0,0)
    S   <- matrix(c(1, tau01, tau01, 1), nrow=2)
    tau <- c(tau0, tau1)
    S   <- diag(tau) %*% S %*% diag(tau)
    U   <- mvrnorm(n, mu=mu, Sigma=S)
    ### simulate AR(1) errors and then the actual outcomes
    dat$eij <- unlist(sapply(p, function(x) arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^2) * sigma))
    dat$logimem_new <- (beta0 + rep(U[,1], times=p)) + (beta1 + rep(U[,2], times=p)) * dat$timesincefirst * dat$group + dat$eij
    dat<-data.frame(dat)
    
    ### subset to prespecified parameters - n, number of visits
    dat <- dat %>% group_by(naccid) %>% filter(n()==minnumvis) %>% ungroup()
    nacc.id<-dat[!duplicated(dat[c("naccid")]),"naccid"]
    id<-data.frame(nacc.id)
    id<-data.frame(sample(id$naccid, sampsize)) 
    colnames(id)[1] <- "naccid"
    N=nrow(id)
    dat<-merge(id,dat, by = "naccid") 
    
    dat.sort <- dat %>% arrange(naccid,vnumber)
    dat.v1 <- dat%>%filter(vnumber==1) %>% mutate(y.mcar = logimem_new) 
    dat.sort2 <- dat.sort  %>% group_by(naccid) %>% filter(row_number()!=1)   %>% ungroup()  
    
    # MCAR
    prop.m = pctmiss  # % missingness
    mcar   = runif(N*(minnumvis-1), min=0, max=1)
    dat.sort2$y.mcar = ifelse(mcar<=prop.m, NA, dat.sort2$logimem_new)  # unrelated to anything
    dat.sort2 <- dat.sort2 %>% arrange(naccid,vnumber)
    dat.miss <- rbind(dat.v1,dat.sort2)
    dat.miss <- dat.miss %>% arrange(naccid,vnumber) 
    nrow(dat.miss)
    dat.miss.check <- dat.miss %>% drop_na("y.mcar")
    nrow(dat.miss.check)
    dat.miss.check <- dat.miss.check  %>% group_by(naccid) %>% filter(n()>2)   %>% slice(1) 
    id.only<- dat.miss.check[,1]
    nrow(dat.miss.check)
    id.only<-as.data.frame(id.only)

    
    dat.misss<-merge(id.only,dat.miss,by="naccid")
    nrow(dat.misss)
    dat.misss <- dat.misss %>% arrange(naccid,vnumber)
    
    
    ### Linear Mixed-Effects, Random Intercept & Random Slope with AR(1) Structure Covariance
    modRIRS_ar1<-lme(y.mcar ~ timesincefirst*group + timesincefirst*baseage , random=~1+timesincefirst|naccid, data=dat.misss, 
                     na.action = na.omit,  correlation = corAR1() ,control = lmeControl(opt = "optim"), method="ML")
    b.lme.rirs.ar1[s]<-summary(modRIRS_ar1)$tTable[5,1]
    se.rirs.ar1[s]<-summary(modRIRS_ar1)$tTable[5,2]
    t05 <- qt(0.975,N-6)
    ci.low.rirs.ar1[s]<-b.lme.rirs.ar1[s] - (t05*se.rirs.ar1[s])
    ci.upp.rirs.ar1[s]<-b.lme.rirs.ar1[s] + (t05*se.rirs.ar1[s])
    tstat.rirs.ar1[s] <- summary(modRIRS_ar1)$tTable[5,4]
    coverage.rirs.ar1[s]<- ((ci.low.rirs.ar1[s] <= truebeta) & (ci.upp.rirs.ar1[s] >= truebeta))
    
    
    ### Linear Mixed-Effects, Random Intercept with AR(1) Structure Covariance
    modRI_ar1<-lme(y.mcar ~ timesincefirst*group +timesincefirst*baseage, random=~1|naccid, data=dat.misss, 
                   na.action = na.omit, correlation = corAR1() ,control = lmeControl(opt = "optim"), method="ML")
    b.lme.ri.ar1[s]<-summary(modRI_ar1)$tTable[5,1]
    se.ri.ar1[s]<-summary(modRI_ar1)$tTable[5,2]
    ci.low.ri.ar1[s]<-b.lme.ri.ar1[s] - (t05*se.ri.ar1[s])
    ci.upp.ri.ar1[s]<-b.lme.ri.ar1[s] + (t05*se.ri.ar1[s])
    tstat.ri.ar1[s] <- summary(modRI_ar1)$tTable[5,4]
    coverage.ri.ar1[s]<- ((ci.low.rirs.ar1[s] <= truebeta) & (ci.upp.rirs.ar1[s] >= truebeta))
    
    ### Linear Mixed-Effects, Random Intercept & Random Slope with Unstructured Covariance
    modRIRS_UN<-lme(y.mcar ~ timesincefirst*group +timesincefirst*baseage, random=~1+timesincefirst|naccid, data=dat.misss ,na.action = na.omit,control = lmeControl(opt = "optim"), method="ML")
    b.lme.rirs.un[s]<-summary(modRIRS_UN)$tTable[5,1]
    se.rirs.un[s]<-summary(modRIRS_UN)$tTable[5,2]
    ci.low.rirs.un[s]<-b.lme.rirs.un[s] - (t05*se.rirs.un[s])
    ci.upp.rirs.un[s]<-b.lme.rirs.un[s] + (t05*se.rirs.un[s])
    tstat.rirs.un[s] <- summary(modRIRS_UN)$tTable[5,4]
    coverage.rirs.un[s]<- ((ci.low.rirs.un[s] <= truebeta) & (ci.upp.rirs.un[s] >= truebeta))
    
    ### Linear Mixed-Effects, Random Intercept with Unstructured Covariance
    modRI_UN<-lme(y.mcar ~ timesincefirst*group + timesincefirst*baseage, random=~1|naccid, data=dat.misss, na.action = na.omit,control = lmeControl(opt = "optim"), method="ML")
    b.lme.ri.un[s]<-summary(modRI_UN)$tTable[5,1]
    se.ri.un[s]<-summary(modRI_UN)$tTable[5,2]
    ci.low.ri.un[s]<-b.lme.ri.un[s] - (t05*se.ri.un[s])
    ci.upp.ri.un[s]<-b.lme.ri.un[s] + (t05*se.ri.un[s])
    tstat.ri.un[s] <- summary(modRI_UN)$tTable[5,4]
    coverage.ri.un[s]<- ((ci.low.ri.un[s] <= truebeta) & (ci.upp.ri.un[s] >= truebeta))
    
    dat.miss2 <- dat.misss %>% arrange(naccid, vnumber) %>% mutate(indx = ifelse(lag(is.na(y.mcar)),1,0)) %>% drop_na(y.mcar) %>%
      mutate(tsl_adj2 = ifelse(indx==1, (timesincefirst - lag(timesincefirst)), tsl_adj))
    dat.miss2 <- dat.miss2 %>% arrange(naccid, vnumber) %>% group_by(naccid) %>% mutate(tsl_a = ifelse(row_number()==1,0,tsl_adj2))
    
    
    dat.age <- reshape2::dcast(dat.miss2, vnumber  ~ naccid, value.var="baseage")
    dat.age <- data.frame(dat.age[,-1])   
    
    # Group 
    dat.grp <- reshape2::dcast(dat.miss2, vnumber  ~ naccid, value.var="group")
    dat.grp <- data.frame(dat.grp[,-1])
    
    #TSL
    dat.tsl <- reshape2::dcast(dat.miss2, vnumber  ~ naccid, value.var="tsl_a")
    dat.tsl <- data.frame(dat.tsl[,-1])
    
    # Logimem Scores
    data_wide <- reshape2::dcast(dat.miss2, vnumber  ~ naccid, value.var="y.mcar")
    data_wide <- data.frame(data_wide[,-1])
    finaldat <- data.frame(cbind(data_wide,dat.grp))
    
    ############## Linear Mixed-Effects Model : RI ###########
    # likes<-rep(0,ncol(data_wide))
    alpha0 <- matrix(c(0,summary(modRIRS_ar1)$tTable[5,1],summary(modRIRS_ar1)$tTable[6,1]), nrow=3, ncol=1 ) # assume alpha_0 = N(0,10^7)
    #sigma0 <- diag(1e7,3)
    sigma0 <- diag(c(1e7,(summary(modRIRS_ar1)$tTable[5,2]^2)*N ,(summary(modRIRS_ar1)$tTable[6,2]^2)*N )) # assume alpha_0 = N(0,10^7)
    #alpha0 <- matrix(c(0,0,0),nrow=3,ncol=1)
    
    tot.like <- 0
    
    fn<-function(par){
      
      for (k in 1:ncol(data_wide)){
        
        age=as.numeric(dat.age[1,k])
        grp=as.numeric(finaldat[1,(k+ncol(data_wide))])
        Gt <- matrix(c(1,0,0,
                       grp,1,0,
                       age,0,1),nrow=3,ncol=3) # transition matrix
        Ft <- matrix(c(1,0,0), nrow=1, ncol=3) # matrix in measurement equation
        W <- diag(c(par[2]^2,0,0))#as.matrix((par[2])) # state variance-covariance matrix
        V <- par[1]^2#as.matrix((par[1]))# measurement variance-covariance matrix
        
        
        yt<-as.matrix(as.numeric(data_wide[,k]))
        yt<-as.matrix(na.omit(yt))
        num<-length(yt)
        #plot.ts(yt)
        
        Gt <- as.matrix(Gt)
        pdim <- nrow(Gt)
        yt <- as.matrix(yt)
        qdim <- ncol(yt)
        
        alphap <- array(NA, dim = c(pdim, 1, num))  # alpha_p= a_{t+1}          
        Pp <- array(NA, dim = c(pdim, pdim, num))  # Pp=P_{t+1}
        alphaf <- array(NA, dim = c(pdim, 1, num))  # alpha_f=a_t
        Pf <- array(NA, dim = c(pdim, pdim, num))  # Pf=P_t
        innov <- array(NA, dim = c(qdim, 1, num))  # innovations
        sig <- array(NA, dim = c(qdim, qdim, num))  # innov var-cov matrix
        Kmat <- array(NA, dim = c(pdim, qdim, num))  # store Kalman gain                  
        
        ########## start filter iterations ###################
        
        # t=1 
        alpha00 <- as.matrix(alpha0, nrow = pdim, ncol = 1)  # state: mean starting value
        P00 <- as.matrix(sigma0, nrow = pdim, ncol = pdim)  # state: variance start value
        
        alphap[, , 1] <- Gt %*% alpha00  # predicted value a_{t+1}
        Pp[, , 1] <- Gt %*% P00 %*% t(Gt) + W  # variance for predicted state value
        sigtemp <- Ft %*% Pp[, , 1] %*% t(Ft) + V  # variance for measurement value
        sig[, , 1] <- (t(sigtemp) + sigtemp)/2  # innov var - make sure it's symmetric
        siginv <- solve(sig[, , 1])  # 1 / innov var 
        K <- Pp[, , 1] %*% t(Ft) %*% siginv  # K_t = predicted variance / innov variance 
        Kmat[, , 1] <- K  # store Kalman gain          
        innov[, , 1] <- yt[1, ] - Ft %*% alphap[, , 1]  # epsilon_t = y_t - a_t
        alphaf[, , 1] <- alphap[, , 1] + K %*% innov[, , 1]  # a_{t+1} = a_t + K_t(epsilon_t)
        Pf[, , 1] <- Pp[, , 1] - Pp[, , 1] %*% t(Ft) %*% t(K) - K %*% Ft %*% Pp[, , 1] + 
          K %*% Ft %*% Pp[, , 1] %*% t(Ft)  %*% t(K) + K %*% sig[, , 1] %*% t(K)
        #Pf[, , 1] <- Pp[, , 1] - K %*% Ft %*% Pp[, , 1]  # variance of forecast
        sigmat <- as.matrix(sig[, , 1], nrow = qdim, ncol = qdim)  # collect variance of measurement errors 
        eigen<-eigen(sigmat)
        cholStatus <- try(u <- chol(sigmat), silent = FALSE)
        cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
        
        # fix the correl matrix
        newMat <- sigmat
        
        
        iter <- 0
        while (cholError) {
          
          iter <- iter + 1
          cat("iteration ", iter, "\n")
          
          # replace -ve eigen values with small +ve number
          newEig <- eigen(newMat)
          newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)
          newEig2 <- as.matrix(newEig2)
          
          # create modified matrix eqn 5 from Brissette et al 2007, inv = transp for
          # eig vectors
          newMat <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)
          
          # normalize modified matrix eqn 6 from Brissette et al 2007
          newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
          
          # try chol again
          cholStatus <- try(u <- chol(newMat), silent = TRUE)
          cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
          
          sigmat<-newMat
        }
        
        
        
        like <-  log(det(sigmat)) + t(innov[, , 1]) %*% siginv %*%
          innov[, , 1]  # calculate -log(likelihood)
        #(num*nrow(Gt)/2)*(log(2*pi)) +
        # t>=2
        
        
        for (i in 2:num) {
          
          vec.tsl <- na.omit(dat.tsl[ ,k])
          
          tsl <-as.numeric(vec.tsl[i])
          
          Gt <- matrix(c(1,0,0,
                         grp*tsl,1,0,
                         age*tsl,0,1),nrow=3,ncol=3) # transition matrix
          
          alphap[, , i] <- Gt %*% alphaf[, , i- 1]  # predicted value a_{t+2}
          Pp[, , i] <- Gt %*% Pf[, , i - 1] %*% t(Gt) + (W*tsl^2)  # variance of predicted state estimate
          
          sigtemp <- Ft %*% Pp[, , i] %*% t(Ft) + V  # variance of measurement error
          sig[, , i] <- (t(sigtemp) + sigtemp)/2  # innov var - make sure it's symmetric
          siginv <- solve(sig[, , i])
          K <- Pp[, , i] %*% t(Ft) %*% siginv  # K_t = predicted variance / innov variance
          Kmat[, , i] <- K  # store Kalman gain          
          innov[, , i] <- yt[i, ] - Ft %*% alphap[, , i]  # epsilon_t = y_t - a_t
          alphaf[, , i] <- alphap[, , i] + K %*% innov[, , i]  # a_{t+1} = a_t + K_t(epsilon_t)
          Pf[, , i] <- Pp[, , i] - Pp[, , i] %*% t(Ft) %*% t(K) - K %*% Ft %*% Pp[, , i] + 
            K %*% Ft %*% Pp[, , i] %*% t(Ft)  %*% t(K) + K %*% sig[, , i] %*% t(K)
          #Pf[, , i] <- Pp[, , i] - K %*% Ft %*% Pp[, , i]  # variance of forecast
          sigmat <- as.matrix(sig[, , i], nrow = qdim, ncol = qdim)  # collect variance of measurement errors
          
          
        }
        like <- like + ((log(det(sigmat)) + t(innov[, , i]) %*% 
                           siginv %*% innov[, , i]))  # calculate -log(likelihood)
        
        #  likes[k] = like
        tot.like <- tot.like + like
      }
      tot.like <- 0.5 * tot.like
      
      return(tot.like)
    }
    
    mle.parms<-optim(c(0.1,0.1), fn,method='L-BFGS-B', hessian=T, control=list(trace=1, REPORT=1))
    mle.parms$convergence # convergence=0 OK!
    
    # Insert Maximized V and W
    var<-rep(0,ncol(data_wide))
    beta.filtered<-rep(0,ncol(data_wide))
    beta.smoothed<-rep(0,ncol(data_wide))
    #sd.beta<-rep(0,nsamp)
    #ci.low<-rep(0,nsamp)
    #ci.upp<-rep(0,nsamp)
    #coverage<-rep(0,nsamp)
    
    for (m in 1:ncol(data_wide)){
      
      #assume alpha_0 = N(0,10^7)# assume alpha_0 = N(0,10^7)
      age=as.numeric(dat.age[1,m])
      grp=as.numeric(finaldat[1,(m+ncol(data_wide))])
      Gt <- matrix(c(1,0,0,
                     grp,1,0,
                     age,0,1),nrow=3,ncol=3) # transition matrix
      Ft <- matrix(c(1,0,0), nrow=1, ncol=3) # matrix in measurement equation
      W <- diag(c(mle.parms$par[2]^2,0,0))#as.matrix((par[2])) # state variance-covariance matrix
      V <- mle.parms$par[1]^2#as.matrix((par[1]))# measurement variance-covariance matrix
      
      
      yt<-as.matrix(as.numeric(data_wide[,m]))
      yt<-as.vector(na.omit(yt))
      num<-length(yt)
      
      #plot.ts(yt)
      
      Gt <- as.matrix(Gt)
      pdim <- nrow(Gt)
      yt <- as.matrix(yt)
      qdim <- ncol(yt)
      
      alphap <- array(NA, dim = c(pdim, 1, num))  # alpha_p= a_{t+1}          
      Pp <- array(NA, dim = c(pdim, pdim, num))  # Pp=P_{t+1}
      alphaf <- array(NA, dim = c(pdim, 1, num))  # alpha_f=a_t
      Pf <- array(NA, dim = c(pdim, pdim, num))  # Pf=P_t
      innov <- array(NA, dim = c(qdim, 1, num))  # innovations
      sig <- array(NA, dim = c(qdim, qdim, num))  # innov var-cov matrix
      Kmat <- array(NA, dim = c(pdim, qdim, num))  # store Kalman gain          
      
      ########## start filter iterations ###################
      
      # t=1 
      alpha00 <- as.matrix(alpha0, nrow = pdim, ncol = 1)  # state: mean starting value
      P00 <- as.matrix(sigma0, nrow = pdim, ncol = pdim)  # state: variance start value
      
      alphap[, , 1] <- Gt %*% alpha00  # predicted value a_{t+1}
      Pp[, , 1] <- Gt %*% P00 %*% t(Gt) + W  # variance for predicted state value
      sigtemp <- Ft %*% Pp[, , 1] %*% t(Ft) + V  # variance for measurement value
      sig[, , 1] <- (t(sigtemp) + sigtemp)/2  # innov var - make sure it's symmetric
      siginv <- solve(sig[, , 1])  # 1 / innov var 
      K <- Pp[, , 1] %*% t(Ft) %*% siginv  # K_t = predicted variance / innov variance 
      Kmat[,,1] <- K  # store Kalman gain          
      innov[, , 1] <- yt[1, ] - Ft %*% alphap[, , 1]  # epsilon_t = y_t - a_t
      alphaf[, , 1] <- alphap[, , 1] + K %*% innov[, , 1]  # a_{t+1} = a_t + K_t(epsilon_t)
      Pf[, , 1] <- Pp[, , 1] - Pp[, , 1] %*% t(Ft) %*% t(K) - K %*% Ft %*% Pp[, , 1] + 
        K %*% Ft %*% Pp[, , 1] %*% t(Ft)  %*% t(K) + K %*% sig[, , 1] %*% t(K)
      #Pf[, , 1] <- Pp[, , 1] - K %*% Ft %*% Pp[, , 1]  # variance of forecast
      sigmat <- as.matrix(sig[, , 1], nrow = qdim, ncol = qdim)  # collect variance of measurement errors 
      
      
      like <- (log(det(sigmat)) + t(innov[, , 1]) %*% siginv %*%
                 innov[, , 1]) # calculate -log(likelihood)
      # t>=2 
      
      
      for (i in 2:num) {
        if (num < 2) 
          break
        
        vec.tsl <- na.omit(dat.tsl[ ,m])
        
        tsl <-as.numeric(vec.tsl[i])
        
        Gt <- matrix(c(1,0,0,
                       grp*tsl,1,0,
                       age*tsl,0,1),nrow=3,ncol=3) # transition matrix
        
        
        alphap[, , i] <- Gt %*% alphaf[, , i- 1] # predicted value a_{t+2}
        Pp[, , i] <- Gt %*% Pf[, , i - 1] %*% t(Gt) + (W*tsl^2)  # variance of predicted state estimate
        
        sigtemp <- Ft %*% Pp[, , i] %*% t(Ft) + V  # variance of measurement error
        sig[, , i] <- (t(sigtemp) + sigtemp)/2  # innov var - make sure it's symmetric
        siginv <- solve(sig[, , i])
        K <- Pp[, , i] %*% t(Ft) %*% siginv  # K_t = predicted variance / innov variance
        Kmat[,,i] <- K  # store Kalman gain          
        innov[, , i] <- yt[i, ] - Ft %*% alphap[, , i]  # epsilon_t = y_t - a_t
        alphaf[, , i] <- alphap[, , i] + K %*% innov[, , i]  # a_{t+1} = a_t + K_t(epsilon_t)
        Pf[, , i] <- Pp[, , i] - Pp[, , i] %*% t(Ft) %*% t(K) - K %*% Ft %*% Pp[, , i] + 
          K %*% Ft %*% Pp[, , i] %*% t(Ft)  %*% t(K) + K %*% sig[, , i] %*% t(K)
        #Pf[, , i] <- Pp[, , i] - K %*% Ft %*% Pp[, , i]  # variance of forecast       
        
        
        kf <- list(alphap = alphap, Pp = Pp, alphaf = alphaf, Pf = Pf, 
                   innov = innov, sig = sig, Kn = K)
        
        ## Kalman smoother
        pdim <- nrow(as.matrix(Gt))
        alphas <- array(NA, dim = c(pdim, 1, num))  # alpha_s = smoothed alpha values 
        Ps <- array(NA, dim = c(pdim, pdim, num))  # Ps=P_t^s
        J <- array(NA, dim = c(pdim, pdim, num))  # J=J_t
        alphas[, , num] <- kf$alphaf[, , num]  # starting value for a_T^s
        Ps[, , num] <- kf$Pf[, , num]  # starting value for prediction variance  
        
        
        ########## start smoother iterations ###################
        w.numert <- rep(0,num) 
        
        for (k in num:2) {
          J[, , k - 1] <- (kf$Pf[, , k - 1] %*% t(Gt)) %*% solve(kf$Pp[, , k])
          alphas[, , k - 1] <- kf$alphaf[, , k - 1] + J[, , k -1] %*% (alphas[, , k] - kf$alphap[, , k])
          Ps[, , k-1] <- kf$Pf[, , k - 1] + J[, , k - 1] %*% (Ps[, , k] - kf$Pp[, , k]) %*% t(J[, , k - 1])
        }
        
        alpha00 <- alpha0
        P00 <- sigma0
        J0 <- as.matrix((P00 %*% t(Gt)) %*% solve(kf$Pp[, , 1]), 
                        nrow = pdim, ncol = pdim)
        alpha0n <- as.matrix(alpha00 + J0 %*% (alphas[, , 1] - kf$alphap[, , 1]), nrow = pdim, ncol = 1)
        P0n <- P00 + J0 %*% (Ps[, , k] - kf$Pp[, , k]) %*% t(J0)
        
        ks <- list(alphas = alphas, Ps = Ps, alpha0n = alpha0n, 
                   P0n = P0n, J0 = J0, J = J, alphap = kf$alphap, Pp = kf$Pp, 
                   alphaf = kf$alphaf, Pf = kf$Pf, like = kf$like, Kn = kf$K)
        
        
        
      }
      
      beta.filtered[m] <- kf$alphap[2,,nrow(yt)] 
      beta.smoothed[m] <- ks$alphas[2,,nrow(yt)]
      
      var[m]<-ks$Ps[2,2,nrow(yt)]
      #sd[m]<-sqrt(ks$Ps[2,2,nrow(dat.outcome)])
      #ci.low[m]<-beta.m[m]-1.96*sd[m]
      #ci.upp[m]<-beta.m[m]+1.96*sd[m]
      
      #coverage[m]<- sum((ci.low[m] <= 0.5) & (ci.upp[m] >= 0.5))
      
      #sd.final<-mean(sd)
      #beta.m.final<-mean(beta.m)
      #F.inv[k]<-as.numeric(siginv)
      #sd[m]<-sqrt(ks$Ps[2,2,nrow(dat.outcome)])
      #ci.low[m]<-beta.m[m]-1.96*sd[m]
      #ci.upp[m]<-beta.m[m]+1.96*sd[m]
      
      #coverage[m]<- sum((ci.low[m] <= 0.5) & (ci.upp[m] >= 0.5))
      
      
      #sd.final<-mean(sd)
      #F.inv[k]<-as.numeric(siginv)
    }
    
    
    
    # collect beta/sd/se/cis/tstat/coverage per simulation
    ###
    beta.llt[s] <-mean(beta.smoothed)
    sd.beta.llt[s] <-sqrt(sum(var)/(N-1))
    
    #sd.beta.llt[s] <-sd(beta.smoothed)
    se.beta.llt[s] <- sd.beta.llt[s]/sqrt(N)
    #sqrt(length(betas$bgroup[group!=0]))
    
    t05 <- qt(0.975,(N-5))
    ci.low.llt[s]<-beta.llt[s]-(t05*se.beta.llt[s])
    ci.upp.llt[s]<-beta.llt[s]+(t05*se.beta.llt[s])
    
    tstat.llt[s] <- beta.llt[s]/se.beta.llt[s]
    
    coverage.llt[s]<- ((ci.low.llt[s] <= truebeta) & (ci.upp.llt[s] >= truebeta))}, error = function(e) {skip_to_next <<- TRUE},
  warning = function(w) {skip_to_next <<- TRUE})

if(skip_to_next) { next }
  }
  list(b.lme.rirs.ar1=b.lme.rirs.ar1, se.rirs.ar1=se.rirs.ar1, ci.low.rirs.ar1=ci.low.rirs.ar1, ci.upp.rirs.ar1=ci.upp.rirs.ar1,coverage.rirs.ar1=coverage.rirs.ar1, tstat.rirs.ar1=tstat.rirs.ar1,
       b.lme.ri.ar1=b.lme.ri.ar1, se.ri.ar1=se.ri.ar1, ci.low.ri.ar1=ci.low.ri.ar1, ci.upp.ri.ar1=ci.upp.ri.ar1,coverage.ri.ar1=coverage.ri.ar1, tstat.ri.ar1=tstat.ri.ar1,
       b.lme.rirs.un=b.lme.rirs.un, se.rirs.un=se.rirs.un, ci.low.rirs.un=ci.low.rirs.un, ci.upp.rirs.un=ci.upp.rirs.un,coverage.rirs.un=coverage.rirs.un, tstat.rirs.un=tstat.rirs.un,
       b.lme.ri.un=b.lme.ri.un, se.ri.un=se.ri.un, ci.low.ri.un=ci.low.ri.un, ci.upp.ri.un=ci.upp.ri.un,coverage.ri.un=coverage.ri.un, tstat.ri.un=tstat.ri.un,
       beta.llt=beta.llt, sd.beta.llt=sd.beta.llt, se.beta.llt=se.beta.llt, ci.low.llt=ci.low.llt, ci.upp.llt=ci.upp.llt,coverage.llt=coverage.llt, tstat.llt=tstat.llt)
}


sim.results<-simulation(1000, 5, 9, 200, 12183, .1)
write.csv(sim.results, "n1000_9vis_beta5_mcar_10.csv")

sim.results<-simulation(1000, 5, 9, 200, 24, .2)
write.csv(sim.results, "n1000_9vis_beta5_mcar_20.csv")

sim.results<-simulation(1000, 5, 9, 200, 4142, .3)
write.csv(sim.results, "n1000_9vis_beta5_mcar_30.csv")

sim.results<-simulation(1000, 5, 9, 200, 678, .4)
write.csv(sim.results, "n1000_9vis_beta5_mcar_40.csv")



sim.results<-simulation(1000, 5, 5, 200, 121683, .1)
write.csv(sim.results, "n1000_5vis_beta5_mcar_10.csv")

sim.results<-simulation(1000, 5, 5, 200, 224, .2)
write.csv(sim.results, "n1000_5vis_beta5_mcar_20.csv")

sim.results<-simulation(1000, 5, 5, 200, 40142, .3)
write.csv(sim.results, "n1000_5vis_beta5_mcar_30.csv")

sim.results<-simulation(1000, 5, 5, 200, 6078, .4)
write.csv(sim.results, "n1000_5vis_beta5_mcar_40.csv")






sim.results<-simulation(1000, .5, 9, 200, 12183, .1)
write.csv(sim.results, "n1000_9vis_beta05_mcar_10.csv")

sim.results<-simulation(1000, .5, 9, 200, 24, .2)
write.csv(sim.results, "n1000_9vis_beta05_mcar_20.csv")

sim.results<-simulation(1000, .5, 9, 200, 4142, .3)
write.csv(sim.results, "n1000_9vis_beta05_mcar_30.csv")

sim.results<-simulation(1000, .5, 9, 200, 678, .4)
write.csv(sim.results, "n1000_9vis_beta05_mcar_40.csv")



sim.results<-simulation(1000, .5, 5, 200, 121683, .1)
write.csv(sim.results, "n1000_5vis_beta05_mcar_10.csv")

sim.results<-simulation(1000, .5, 5, 200, 224, .2)
write.csv(sim.results, "n1000_5vis_beta05_mcar_20.csv")

sim.results<-simulation(1000, .5, 5, 200, 40142, .3)
write.csv(sim.results, "n1000_5vis_beta05_mcar_30.csv")

sim.results<-simulation(1000, .5, 5, 200, 6078, .4)
write.csv(sim.results, "n1000_5vis_beta05_mcar_40.csv")










sim.results<-simulation(1000, 5, 9, 100, 12183, .1)
write.csv(sim.results, "n1000_9vis_beta5_mcar_10_100.csv")

sim.results<-simulation(1000, 5, 9, 100, 24, .2)
write.csv(sim.results, "n1000_9vis_beta5_mcar_20_100.csv")

sim.results<-simulation(1000, 5, 9, 100, 4142, .3)
write.csv(sim.results, "n1000_9vis_beta5_mcar_30_100.csv")

sim.results<-simulation(1000, 5, 9, 100, 678, .4)
write.csv(sim.results, "n1000_9vis_beta5_mcar_40_100.csv")



sim.results<-simulation(1000, 5, 5, 100, 121683, .1)
write.csv(sim.results, "n1000_5vis_beta5_mcar_10_100.csv")

sim.results<-simulation(1000, 5, 5, 100, 224, .2)
write.csv(sim.results, "n1000_5vis_beta5_mcar_20_100.csv")

sim.results<-simulation(1000, 5, 5, 100, 40142, .3)
write.csv(sim.results, "n1000_5vis_beta5_mcar_30_100.csv")

sim.results<-simulation(1000, 5, 5, 100, 6078, .4)
write.csv(sim.results, "n1000_5vis_beta5_mcar_40_100.csv")







sim.results<-simulation(1000, .5, 9, 100, 12183, .1)
write.csv(sim.results, "n1000_9vis_beta05_mcar_10_100.csv")

sim.results<-simulation(1000, .5, 9, 100, 24, .2)
write.csv(sim.results, "n1000_9vis_beta05_mcar_20_100.csv")

sim.results<-simulation(1000, .5, 9, 100, 4142, .3)
write.csv(sim.results, "n1000_9vis_beta05_mcar_30_100.csv")

sim.results<-simulation(1000, .5, 9, 100, 678, .4)
write.csv(sim.results, "n1000_9vis_beta05_mcar_40_100.csv")



sim.results<-simulation(1000, .5, 5, 100, 121683, .1)
write.csv(sim.results, "n1000_5vis_beta05_mcar_10_100.csv")

sim.results<-simulation(1000, .5, 5, 100, 224, .2)
write.csv(sim.results, "n1000_5vis_beta05_mcar_20_100.csv")

sim.results<-simulation(1000, .5, 5, 100, 40142, .3)
write.csv(sim.results, "n1000_5vis_beta05_mcar_30_100.csv")

sim.results<-simulation(1000, .5, 5, 100, 6078, .4)
write.csv(sim.results, "n1000_5vis_beta05_mcar_40_100.csv")






sim.results<-simulation(1000, 5, 9, 50, 12183, .1)
write.csv(sim.results, "n1000_9vis_beta5_mcar_10_50.csv")

sim.results<-simulation(1000, 5, 9, 50, 24, .2)
write.csv(sim.results, "n1000_9vis_beta5_mcar_20_50.csv")

sim.results<-simulation(1000, 5, 9, 50, 4142, .3)
write.csv(sim.results, "n1000_9vis_beta5_mcar_30_50.csv")

sim.results<-simulation(1000, 5, 9, 50, 678, .4)
write.csv(sim.results, "n1000_9vis_beta5_mcar_40_50.csv")



sim.results<-simulation(1000, 5, 5, 50, 121683, .1)
write.csv(sim.results, "n1000_5vis_beta5_mcar_10_50.csv")

sim.results<-simulation(1000, 5, 5, 50, 224, .2)
write.csv(sim.results, "n1000_5vis_beta5_mcar_20_50.csv")

sim.results<-simulation(1000, 5, 5, 50, 40142, .3)
write.csv(sim.results, "n1000_5vis_beta5_mcar_30_50.csv")

sim.results<-simulation(1000, 5, 5, 50, 6078, .4)
write.csv(sim.results, "n1000_5vis_beta5_mcar_40_50.csv")






sim.results<-simulation(1000, .5, 9, 50, 12183, .1)
write.csv(sim.results, "n1000_9vis_beta05_mcar_10_50.csv")

sim.results<-simulation(1000, .5, 9, 50, 24, .2)
write.csv(sim.results, "n1000_9vis_beta05_mcar_20_50.csv")

sim.results<-simulation(1000, .5, 9, 50, 4142, .3)
write.csv(sim.results, "n1000_9vis_beta05_mcar_30_50.csv")

sim.results<-simulation(1000, .5, 9, 50, 678, .4)
write.csv(sim.results, "n1000_9vis_beta05_mcar_40_50.csv")



sim.results<-simulation(1000, .5, 5, 50, 121683, .1)
write.csv(sim.results, "n1000_5vis_beta05_mcar_10_50.csv")

sim.results<-simulation(1000, .5, 5, 50, 224, .2)
write.csv(sim.results, "n1000_5vis_beta05_mcar_20_50.csv")

sim.results<-simulation(1000, .5, 5, 50, 40142, .3)
write.csv(sim.results, "n1000_5vis_beta05_mcar_30_50.csv")

sim.results<-simulation(1000, .5, 5, 50, 6078, .4)
write.csv(sim.results, "n1000_5vis_beta05_mcar_40_50.csv")



