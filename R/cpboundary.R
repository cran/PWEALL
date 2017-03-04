cpboundary<-function(Dplan=300,alpha=0.05,two.sided=1,pi1=0.5,cpcut=c(0.2,0.3,0.4),BetaD=log(0.8),Beta0=log(1),prop=seq(0.1,0.9,by=0.1)){
       nr<-length(prop);nc<-length(cpcut)
       alpha1<-0.05/2
       pi0<-1-pi1
       if (two.sided!=1) alpha1<-alpha
       CPNbound<-CPTbound<-CPDbound<-matrix(0,nrow=nr,ncol=nc)
       for (i in 1:nr){
           for (j in 1:nc){
               CPTbound[i,j]<-(-qnorm(cpcut[j])*sqrt(1-prop[i])-qnorm(1-alpha1))/sqrt(pi1*pi0*Dplan)
               CPNbound[i,j]<-(-qnorm(cpcut[j])*sqrt(1-prop[i])-qnorm(1-alpha1))/sqrt(pi1*pi0*Dplan)/prop[i]
               CPDbound[i,j]<-((-qnorm(cpcut[j])*sqrt(1-prop[i])-qnorm(1-alpha1))/sqrt(pi1*pi0*Dplan)-BetaD*(1-prop[i]))/prop[i]
           }
       }
       cnames<-rep("CP",nc)
       for (j in 1:nc)cnames[j]<-paste("CP=",cpcut[j])
       eCPTbound<-data.frame(exp(CPTbound))
       eCPNbound<-data.frame(exp(CPNbound))
       eCPDbound<-data.frame(exp(CPDbound))
       colnames(eCPTbound)<-cnames
       colnames(eCPNbound)<-cnames
       colnames(eCPDbound)<-cnames
       list(CPTbound=eCPTbound,CPNbound=eCPNbound,CPDbound=eCPDbound)
       #The boundary means if the obs. HR is above the boundary, the CP will fall below the designated level at cpcut
}
