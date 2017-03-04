cp<-function(Dplan=300,alpha=0.05,two.sided=1,pi1=0.5,Obsbeta=log(seq(1,0.6,by=-0.01)),BetaD=log(0.8),Beta0=log(1),prop=seq(0.1,0.9,by=0.1)){
       nr<-length(prop);nco<-length(Obsbeta)
       alpha1<-0.05/2
       pi0<-1-pi1
       if (two.sided!=1) alpha1<-alpha
       CPN<-CPT<-CPD<-matrix(0,nrow=nr,ncol=nco)
       for (i in 1:nr){
            for (j in 1:nco){
                CPT[i,j]<-pnorm(-qnorm(1-alpha1)/sqrt(1-prop[i])-sqrt(pi1*pi0*Dplan)*Obsbeta[j]/sqrt(1-prop[i]))
                CPN[i,j]<-pnorm(-qnorm(1-alpha1)/sqrt(1-prop[i])-sqrt(pi1*pi0*Dplan)*Obsbeta[j]*prop[i]/sqrt(1-prop[i]))
                CPD[i,j]<-pnorm(-qnorm(1-alpha1)/sqrt(1-prop[i])-sqrt(pi1*pi0*Dplan)*(Obsbeta[j]*prop[i]+BetaD*(1-prop[i]))/sqrt(1-prop[i]))
            }
       }
       cnames<-rep("0",nco)
       for (j in 1:nco)cnames[j]<-paste("HR=",round(exp(Obsbeta[j]),digits=3))
       colnames(CPT)<-cnames
       colnames(CPN)<-cnames
       colnames(CPD)<-cnames
       list(CPT=CPT,CPN=CPN,CPD=CPD)
       #Conditional power given the obs. hazard ratio=exp(Obsbeta)
}





























