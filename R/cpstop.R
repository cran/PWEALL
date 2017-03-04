cpstop<-function(Dplan=300,pi1=0.5,Beta1=log(0.8),Beta0=log(1),prop=seq(0.1,0.9,by=0.1),HRbound=rep(0.85,length(prop))){
       nr<-length(prop);nc<-length(HRbound)
       pi0<-1-pi1
       pstop0<-pstop1<-rep(0,nr)
       for (i in 1:nr){
               pstop0[i]<-1-pnorm((log(HRbound[i])-Beta0)*sqrt(pi1*pi0*prop[i]*Dplan))
               pstop1[i]<-1-pnorm((log(HRbound[i])-Beta1)*sqrt(pi1*pi0*prop[i]*Dplan))
       }
       list(pstop0=pstop0,pstop1=pstop1)
       #Probability of stopping given the HR boundary (HRbound) at each of the proportion of the number of target events
}




