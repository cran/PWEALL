##########################################################################################
#  random number generator for PieceWise Exp dist with CROSSOVER effect
#  version 1.0 (11/04/2016)
#  version 2.0 (10/10/2017)
#  version 3.0 (10/17/2018): change the algorithm, add two more types of hybrid crossover and output crossover indicators
##########################################################################################
#rPWECROSSFRx<-function(nr=1,rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),rate4=rate2,rate5=rate2,tchange=c(0,1),type=1,ISEEDS=c(13,17,19,23,29)){
rpwecx<-function (nr = 1, rate1 = c(1, 0.5), rate2 = rate1, rate3 = c(0.7,0.4), rate4 = rate2, rate5 = rate2, tchange = c(0, 1), type = 1, 
                   rp2 = 0.5){
  # rate1: hazard before crossover
  # rate2: hazard after crossover
  # rate3: hazard for crossover
  # rate4: hazard after crossover that will be combined with rate2
  # rate5: hazard after crossover that will be combined with rate2
  # tchange: the time points at which either rate1-rate5 changes
  # type: type of crossover, i.e. markov, semi-markov,hybird
  # rp2: re-randomization probability to receive the rescue treatment when semi-markov crossover occurs.
  #      When it happens, the overall hazard will be rp2*rate2(t-s)+(1-rp2)*rate4(t), where rate2 is the hazard 
  #           for the (semi-markov) rescue treatment and rate4 is hazard for the (markov) rescue treatment
  
  if (min(rate3) < 1e-08) 
    rate3 <- rate3 + 1e-08
  ###indicators for crossover
  cxind<-matrix(0,nrow=nr,ncol=2)
  
  y3 <- rpwe(nr = nr, rate = rate3, tchange = tchange)$r
  y1 <- rpwe(nr = nr, rate = rate1, tchange = tchange)$r
  cxind[,1]=cxind[,2]=ifelse(y3<y1,1,0)
  if (type == 1) {
    s2y3 <- pwe(t = y3, rate = rate2, tchange = tchange)$surv
    x <- runif(nr)
    x1 <- x * s2y3
    y <- rep(0, nr)
    y[y3>=y1] <- y1[y3>=y1]
    y[y3<y1] <- qpwe(p = 1 - x1[y3<y1],rate = rate2,tchange = tchange)$q
  }
  else if (type == 2) {
    x <- runif(nr)
    y <- rep(0, nr)
    y <- rep(0, nr)
    y[y3>=y1] <- y1[y3>=y1]
    y[y3<y1] <- y3[y3<y1]+qpwe(p = 1 - x[y3<y1],rate = rate2,tchange = tchange)$q
  }
  else if (type == 3) {
    s4y3 <- pwe(t = y3, rate = (1 - rp2) * rate4, tchange = tchange)$surv
    x <- runif(nr)
    z <- runif(nr)
    z1 <- z * s4y3
    y <- rep(0, nr)
    y <- rep(0, nr)
    y[y3>=y1] <- y1[y3>=y1]
    if (sum(y3<y1) > 0) {
      atemp <- y3[y3<y1] + qpwe(p = 1 - x[y3<y1], 
                                rate = rp2 * rate2, tchange = tchange)$q
      btemp <- qpwe(p = 1 - z1[y3<y1], rate = (1 - 
                                                 rp2) * rate4, tchange = tchange)$q
      y[y3<y1] <- pmin(atemp, btemp)
      cxind[y3<y1,2]=cxind[y3<y1,1]*ifelse(atemp<=btemp,1,0)
    }
  }
  else if (type == 4) {
    s4y3 <- pwe(t = y3, rate = rate4, tchange = tchange)$surv
    x <- runif(nr)
    x1 <-x*s4y3
    z <- rbinom(nr,size=1,prob=rp2)
    cxind[,2]=z*cxind[,1]
    y <- rep(0, nr)
    y[y3>=y1] <-y1[y3>=y1]
    y[y3<y1&z==1]<-y3[y3<y1&z==1]+qpwe(p = 1 - x[y3<y1&z==1],rate = rate2,tchange = tchange)$q
    y[y3<y1&z==0]<-qpwe(p = 1 - x1[y3<y1&z==0],rate = rate4,tchange = tchange)$q
  }
  else if (type == 5) {
    s4y3 <- pwe(t = y3, rate = rate4, tchange = tchange)$surv
    x <- runif(nr)
    x1 <-x*s4y3
    z <- rbinom(nr,size=1,prob=rp2)
    cxind[,2]=z*cxind[,1]
    y <- rep(0, nr)
    y[y3>=y1] <-y1[y3>=y1]
    y[y3<y1&z==1]<-y3[y3<y1&z==1]+qpwe(p = 1 - x[y3<y1&z==1],rate = rate2,tchange = tchange)$q
    y[y3<y1&z==0]<-y3[y3<y1&z==0]+qpwe(p = 1 - x[y3<y1&z==0],rate = rate4,tchange = tchange)$q
  }
  list(r=y,rx=y3,cxind=cxind,type=type)
}



