overallcov<-function(tfix=2.0,tfix0=1.0,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
              rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
              rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
              rate10=c(1,0.7),rate20=rate10,rate30=rate31,
              rate40=rate20,rate50=rate20,ratec0=ratec1,
              tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,
              eps=1.0e-2,veps=1.0e-2,beta=0,beta0=0){

  part1<-overallcovp1(tfix=tfix,tfix0=tfix0,taur=taur,u=u,ut=ut,pi1=pi1,
              rate11=rate11,rate21=rate21,rate31=rate31,
              rate41=rate41,rate51=rate51,ratec1=ratec1,
              rate10=rate10,rate20=rate20,rate30=rate30,
              rate40=rate40,rate50=rate50,ratec0=ratec0,
              tchange=tchange,type1=type1,type0=type0,rp21=rp21,rp20=rp20,
              eps=eps,veps=veps,beta=beta,beta0=beta0)

  part2<-overallcovp2(tfix=tfix,tfix0=tfix0,taur=taur,u=u,ut=ut,pi1=pi1,
              rate11=rate11,rate21=rate21,rate31=rate31,
              rate41=rate41,rate51=rate51,ratec1=ratec1,
              rate10=rate10,rate20=rate20,rate30=rate30,
              rate40=rate40,rate50=rate50,ratec0=ratec0,
              tchange=tchange,type1=type1,type0=type0,rp21=rp21,rp20=rp20,
              eps=eps,veps=veps,beta=beta,beta0=beta0)

  EA1<-part1$EA1
  EA2<-part2$EA2
  covbeta<-part1$covbeta1+part2$cov234
  covbeta1<-part1$covbeta1
  covbeta2<-part2$covbeta2
  covbeta3<-part2$covbeta3
  covbeta4<-part2$covbeta4
  list(covbeta=covbeta,covbeta1=covbeta1,covbeta2=covbeta2,covbeta3=covbeta3,covbeta4=covbeta4,EA1=EA1,EA2=EA2)
}
