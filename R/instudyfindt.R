#################################################################################
##   This function will provide in-study timeline prediction 
##   version 1.0 (1/19/2017)
#################################################################################

instudyfindt<-function(target=400,y=exp(rnorm(300)),z=rbinom(300,1,0.5),d=rep(c(0,1,2),each=100),
                        tcut=2,blinded=1,type0=1,type1=type0,rp20=0.5,rp21=0.5,tchange=c(0,1),
                        rate10=c(1,0.7),rate20=c(0.9,0.7),rate30=c(0.4,0.6),rate40=rate20,
                        rate50=rate20,ratec0=c(0.3,0.3),
                        rate11=rate10,rate21=rate20,rate31=rate30,
                        rate41=rate40,rate51=rate50,ratec1=ratec0,
                        withmorerec=1,
                        ntotal=1000,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                        ntype0=1,ntype1=1,nrp20=0.5,nrp21=0.5,ntchange=c(0,1),
                        nrate10=rate10,nrate20=rate20,nrate30=rate30,nrate40=rate40,
                        nrate50=rate50,nratec0=ratec0,
                        nrate11=rate10,nrate21=rate20,nrate31=rate30,nrate41=rate40,
                        nrate51=rate50,nratec1=ratec0,
                        eps=1.0e-2,init=tcut*1.1,epsilon=0.001,maxiter=100){
  ieps<-1
  iter<-0
  ny<-length(y)
  t0<-init;nd1<-sum(d==1);nd2<-sum(d==2)
  dvprimehist<-dvaluehist<-t1hist<-rep(0,maxiter);
  while (ieps>epsilon & iter<=maxiter){
    iter<-iter+1
    dvalue<-nd1;d2value<-dvprime<-0

    ##Handling the exisiting data
    if (nd2>0){
      y2<-y[d==2]
      z2<-z[d==2]
      n21<-sum(z2==1)
      n20<-sum(z2==0)
      at1<-(t0-tcut)+y2
      at0<-y2
      if (blinded==1){
            a1<-pwecxcens(t=c(at1),rate1=rate10,rate2=rate20,rate3=rate30,rate4=rate40,rate5=rate50,ratec=ratec0,
                          tchange=tchange,type=type0,rp2=rp20,eps=eps)
            a0<-pwecxcens(t=c(at0),rate1=rate10,rate2=rate20,rate3=rate30,rate4=rate40,rate5=rate50,ratec=ratec0,
                    tchange=tchange,type=type0,rp2=rp20,eps=eps)
            ap<-(a1$du-a0$du)/(a0$s*a0$sc)
            dvalue<-dvalue+sum(ap)
            dvprime<-dvprime+sum(a1$duprime/(a0$s*a0$sc))
            d2value<-d2value+sum(ap*(1-ap))
      }
      else if (blinded==0){
          if (n20>0){
              a1<-pwecxcens(t=c(at1[z2==0]),rate1=rate10,rate2=rate20,rate3=rate30,rate4=rate40,rate5=rate50,ratec=ratec0,
                      tchange=tchange,type=type0,rp2=rp20,eps=eps)
              a0<-pwecxcens(t=c(at0[z2==0]),rate1=rate10,rate2=rate20,rate3=rate30,rate4=rate40,rate5=rate50,ratec=ratec0,
                      tchange=tchange,type=type0,rp2=rp20,eps=eps)
              ap<-(a1$du-a0$du)/(a0$s*a0$sc)
              dvalue<-dvalue+sum(ap)
              dvprime<-dvprime+sum(a1$duprime/(a0$s*a0$sc))
              d2value<-d2value+sum(ap*(1-ap))
          }
          if (n21>0){
              a1<-pwecxcens(t=c(at1[z2==1]),rate1=rate11,rate2=rate21,rate3=rate31,rate4=rate41,rate5=rate51,ratec=ratec1,
                        tchange=tchange,type=type1,rp2=rp21,eps=eps)
              a0<-pwecxcens(t=c(at0[z2==1]),rate1=rate11,rate2=rate21,rate3=rate31,rate4=rate41,rate5=rate51,ratec=ratec1,
                        tchange=tchange,type=type1,rp2=rp21,eps=eps)
              ap<-(a1$du-a0$du)/(a0$s*a0$sc)
              dvalue<-dvalue+sum(ap)
              dvprime<-dvprime+sum(a1$duprime/(a0$s*a0$sc))
              d2value<-d2value+sum(ap*(1-ap))
          }
      }
    }
    
    ##Handling the new data
    if (withmorerec==1){
        dp1<-pwecxpwu(t=c(t0-tcut),taur=taur,u=u,ut=ut,rate1=nrate11,rate2=nrate21,rate3=nrate31,
                         rate4=nrate41,rate5=nrate51,ratec=nratec1,tchange=ntchange,type=ntype1,rp2=nrp21,eps=eps)
        dp0<-pwecxpwu(t=c(t0-tcut),taur=taur,u=u,ut=ut,rate1=nrate10,rate2=nrate20,rate3=nrate30,
                         rate4=nrate40,rate5=nrate50,ratec=nratec0,tchange=ntchange,type=ntype0,rp2=nrp20,eps=eps)

        dvalue<-dvalue+ntotal*pi1*dp1$du+ntotal*(1-pi1)*dp0$du
        dvprime<-dvprime+ntotal*pi1*dp1$duprime+ntotal*(1-pi1)*dp0$duprime
        ny<-ny+ntotal
        d2value<-d2value+ntotal*pi1*dp1$du*(1-dp1$du)+ntotal*(1-pi1)*dp0$du*(1-dp0$du)
    }

    ###Calculate the next step $t1$
    dvprime1<-dvprime
    enn<-0.0001*ny
    if (dvprime>=0.0 & dvprime<=enn)dvprime1<-enn
    else if (dvprime<0 & dvprime>=-enn) dvprime1<--enn
             
    t1<-t0-(dvalue-target)/dvprime1
    ieps<-abs(t1-t0)
    t0<-t1
    t1hist[iter]<-t1
    dvaluehist[iter]<-dvalue
    dvprimehist[iter]<-dvprime
    tvar<-ny*d2value/dvprime^2
  }
  list(t1=t1,dvalue=dvalue,dvprime=dvprime,tvar=tvar,ny=ny,eps=ieps,iter=iter,t1hist=t1hist[1:iter],dvaluehist=dvaluehist[1:iter],dvprimehist=dvprimehist[1:iter])
}
