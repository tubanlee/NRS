#because I combined all estimators into one function, there might be errors when running several functions in a time and these are not bugs,
#but because R is prone to produce errors for such a large function (I haven't found anything wrong), 
#Run one function one time. If there is an error, try to restart and then it will be fixed.
#It will be completely fixed in the future by writing the codes in C++.
library(Lmoments)
least_common_multiple<-function(a,b){
  g<-greatest_common_divisor(a, b)
  (a/g * b)
}
greatest_common_divisor<- function(a, b) {
  if (b == 0) a else Recall(b, a %% b)
}
data_augmentation<-function (x,targetsize){
  lengthx<-length(x)
  meanx<-mean(x)
  disn<-least_common_multiple(lengthx,targetsize)
  orderedx<-sort(x,decreasing = FALSE,method ="radix")
  disori<-rep(orderedx, each = disn/lengthx)
  group<-rep(1:targetsize, each =disn/targetsize)
  data_augmentationresult<-sapply(split(disori, group), mean)
  data_augmentationresult
}
#equinveral trimmed mean and complementary trimmed mean
etm<-function (x,interval=9,fast=TRUE,batch=1000){
  lengthx<-length(x)
  Ksamples<-lengthx/interval
  IntKsamples<-ceiling(Ksamples)
  target1<-IntKsamples*interval
  if (Ksamples%%1!=0 ){
    if (fast){
      if(lengthx<10000){
        x_ordered<-data_augmentation(x,target1)
      }
      else{
        x<-c(sample(x,target1-lengthx,replace = FALSE),x)
        x_ordered<-sort(x,decreasing = FALSE,method ="radix")
      }
      group<-rep(rep(c(1,2,3), each=IntKsamples), times=target1/(IntKsamples*3))
      group<-replace(group,(target1-IntKsamples+1):target1, 0)
      group<-replace(group,1:(IntKsamples), 0)
      group[group == 1] <- 3
      Group1<-split(x_ordered, group)
      Groupsum<-(sapply(Group1, sum))
      partiallist<-(c(x_ordered[(IntKsamples+1):(2*(IntKsamples))],x_ordered[(target1-2*IntKsamples+1):(target1-IntKsamples)]))
      Groupsumweight<-c(Groupsum[1],(sum(c(Groupsum[2],sum(partiallist)))),Groupsum[3])
      etmlength<-length(which(group == 2))+length(partiallist)
      ctmlength<-length(which(group == 3))
      Groupmean<-c(Groupsumweight[1]/length(which(group == 0)),Groupsumweight[2]/etmlength,Groupsumweight[3]/ctmlength)
      (Groupmean)
    }
    else{
      batch1<-c()
      batch2<-c()
      batch3<-c()
      for (i in 1:batch){
        x1<-c(sample(x,target1-lengthx,replace = FALSE),x)
        x_ordered<-sort(x1,decreasing = FALSE,method ="radix")
        group<-rep(rep(c(1,2,3), each=IntKsamples), times=target1/(IntKsamples*3))
        group<-replace(group,(target1-IntKsamples+1):target1, 0)
        group<-replace(group,1:(IntKsamples), 0)
        group[group == 1] <- 3
        Group1<-split(x_ordered, group)
        Groupsum<-(sapply(Group1, sum))
        partiallist<-(c(x_ordered[(IntKsamples+1):(2*(IntKsamples))],x_ordered[(target1-2*IntKsamples+1):(target1-IntKsamples)]))
        Groupsumweight<-c(Groupsum[1],(sum(c(Groupsum[2],sum(partiallist)))),Groupsum[3])
        etmlength<-length(which(group == 2))+length(partiallist)
        ctmlength<-length(which(group == 3))
        Groupmean<-c(Groupsumweight[1]/length(which(group == 0)),Groupsumweight[2]/etmlength,Groupsumweight[3]/ctmlength)
        batch1<-c(batch1,Groupmean[1])
        batch2<-c(batch2,Groupmean[2])
        batch3<-c(batch3,Groupmean[3])
      }
      ((c((mean(batch1)),mean(batch2),mean(batch3))))
    }
  }
  else{
    x_ordered<-sort(x,decreasing = FALSE,method ="radix")
    group<-rep(rep(c(1,2,3), each=IntKsamples), times=lengthx/(IntKsamples*3))
    group<-replace(group,(target1-IntKsamples+1):target1, 0)
    group<-replace(group,1:(IntKsamples), 0)
    group[group == 1] <- 3
    Group1<-split(x_ordered, group)
    Groupsum<-(sapply(Group1, sum))
    partiallist<-(c(x_ordered[(IntKsamples+1):(2*(IntKsamples))],x_ordered[(target1-2*IntKsamples+1):(target1-IntKsamples)]))
    Groupsumweight<-c(Groupsum[1],(sum(c(Groupsum[2],sum(partiallist)))),Groupsum[3])
    etmlength<-length(which(group == 2))+length(partiallist)
    ctmlength<-length(which(group == 3))
    Groupmean<-c(Groupsumweight[1]/length(which(group == 0)),Groupsumweight[2]/etmlength,Groupsumweight[3]/ctmlength)
    (Groupmean)}
}

mmm<-function(x,interval=9,fast=TRUE,batch=1000,sorted=FALSE,drm=0.3665,dqm=0.82224){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  etm1<-etm(x,interval=interval,fast=fast,batch=batch)
  mx1<-(min(which(sortedx>(etm1[2])))-1)/length(x)
  mx2<-1/2
  if (mx1>0.5){
    quatiletarget<-abs(1-mx1)*((mx1-mx2)*2)*(((abs(mx1-mx2)*2))^dqm)+mx1
  }else{
    quatiletarget<-abs(0-mx1)*((mx1-mx2)*2)*(((abs(mx1-mx2)*2))^dqm)+mx1
  }
  qm1<-quantile(sortedx,quatiletarget)
  
  rm1<--drm*etm1[3]+etm1[2]+drm*etm1[2]
  listd<-c(rm1,qm1)
  return(listd)
}

mmme<-function(x,interval=9,fast=TRUE,batch=1000,sorted=FALSE,drm=0.3665,dqm=0.82224){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  etm1<-etm(x,interval=interval,fast=fast,batch=batch)
  mx1<-(min(which(sortedx>(etm1[2])))-1)/length(x)
  mx2<-1/2
  if (mx1>0.5){
    quatiletarget<-abs(1-mx1)*((mx1-mx2)*2)*(((abs(mx1-mx2)*2))^dqm)+mx1
  }else{
    quatiletarget<-abs(0-mx1)*((mx1-mx2)*2)*(((abs(mx1-mx2)*2))^dqm)+mx1
  }
  qm1<-quantile(sortedx,quatiletarget)
  rm1<--drm*etm1[3]+etm1[2]+drm*etm1[2]
  listd<-c(mean1=mean(sortedx),etm1=etm1[2],rm1=rm1,qm1=qm1)
  return(listd)
}

rqscale<-function (x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,dlrm=0.3659,drm=0.7930,dlqm=0.8218,dqm=0.7825){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  if (boot){
    subtract<-t(replicate(subsample, sort(sample(sortedx, size = 2))))
    getlm<-function(vector){ 
      (vector[2]-vector[1])/2
    }
    alllm<-function(sortedx){ 
      subtract<-t(combn(sortedx, 2))
      apply(subtract,MARGIN=1,FUN=getlm)
    }
    dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
    getm<-function(vector){ 
      ((vector[1]-vector[2])^2)/2
    }
    allm<-function(sortedx){ 
      subtract<-t(combn(sortedx, 2))
      apply(subtract,MARGIN=1,FUN=getm)
    }
    dp2m<-apply(subtract,MARGIN=1,FUN=getm)
  }else{
    subtract<-t(combn(sortedx, 2))
    subtract<-sapply(sortedx, "-", sortedx)
    subtract[lower.tri(subtract)] <- NA
    diag(subtract)=NA
    subtract<-na.omit(as.vector(subtract))
    dp<-subtract[subtract>0]
    dp2lm<-dp/2
    dp2m<-(dp^2)/2
  }
  lengthn<-length(sortedx)
  lm1<-Lmoments(sortedx)
  expectdps<-lm1[2]
  expectdp2s<-(sd(sortedx))^2
  dlmo<-mmm(x=dp2lm,interval=9,fast=TRUE,batch=10000,sorted=FALSE,drm=dlrm,dqm=dlqm)
  dmo<-mmm(x=dp2m,interval=9,fast=TRUE,batch=10000,sorted=FALSE,drm=drm,dqm=dqm)
  all<-c(l2=expectdps,rl2=dlmo[1],ql2=dlmo[2],sdl2=sd(dp2lm),
         var=expectdp2s,rvar=dmo[1],qvar=dmo[2],sdvar=sd(dp2m)
  )
  return(all)
}


rqtm<-function (x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,dlrm=0.1808,drm=1.7492,dlqm=1.1753,dqm=0.5715){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  if (boot){
    subtract<-t(replicate(subsample, sort(sample(sortedx, size = 3))))
  }else{
    subtract<-t(combn(sortedx, 3))
  }
  getlm<-function(vector){ 
    (1/3)*(vector[3]-2*vector[2]+vector[1])
  }
  alllm<-function(sortedx){ 
    subtract<-t(combn(sortedx, 3))
    apply(subtract,MARGIN=1,FUN=getlm)
  }
  dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
  
  getm<-function(vector){ 
    ((1/6)*(2*vector[1]-vector[2]-vector[3])*(-1*vector[1]+2*vector[2]-vector[3])*(-vector[1]-vector[2]+2*vector[3]))
  }
  allm<-function(sortedx){ 
    subtract<-t(combn(sortedx, 3))
    apply(subtract,MARGIN=1,FUN=getm)
  }
  dp2m<-apply(subtract,MARGIN=1,FUN=getm)
  
  lengthn<-length(sortedx)
  lm1<-Lmoments(sortedx)
  
  expectdps<-lm1[3]
  expectdp2s<-((sum((sortedx - mean(sortedx))^3)/lengthn)*(lengthn^2/((lengthn-1)*(lengthn-2))))
  dlmo<-mmm(x=dp2lm,interval=9,fast=fast,batch=batch,sorted=FALSE,drm=dlrm,dqm=dlqm)
  dmo<-mmm(x=dp2m,interval=9,fast=fast,batch=batch,sorted=FALSE,drm=drm,dqm=dqm)
  
  all<-c(l3=expectdps,rl3=dlmo[1],ql3=dlmo[2],sdl3=sd(dp2lm),
         tm=expectdp2s,rtm=dmo[1],qtm=dmo[2],sdtm=sd(dp2m)
  )
  return(all)
}


rqfm<-function (x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,dlrm=-0.3542,drm=3.4560,dlqm=NaN,dqm=0.1246){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  getm<-function(vector){ 
    resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                  4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                  vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                  6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                  vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                               6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
    (resd)
  }
  allm<-function(sortedx){ 
    subtract<-t(combn(sortedx, 4))
    apply(subtract,MARGIN=1,FUN=getm)
  }
  if (boot){
    subtract<-t(replicate(subsample, sort(sample(sortedx, size = 4))))
  }else{
    subtract<-t(combn(sortedx, 4))
  }
  
  dp2m<-apply(subtract,MARGIN=1,FUN=getm)
  
  getlm<-function(vector){ 
    (1/4)*(vector[4]-3*vector[3]+3*vector[2]-vector[1])
  }
  alllm<-function(sortedx){ 
    subtract<-t(combn(sortedx, 4))
    apply(subtract,MARGIN=1,FUN=getlm)
  }
  
  dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
  
  lengthn<-length(sortedx)
  lm1<-Lmoments(sortedx)
  expectdps<-lm1[4]
  expectdp2s<-(sum((sortedx - mean(sortedx))^4)/lengthn)
  dlmo<-mmm(x=dp2lm,interval=9,fast=fast,batch=batch,sorted=FALSE,drm=dlrm,dqm=dlqm)
  dmo<-mmm(x=dp2m,interval=9,fast=fast,batch=batch,sorted=FALSE,drm=drm,dqm=dqm)

  all<-c(l4=expectdps,rl4=dlmo[1],ql4=dlmo[2],sdl4=sd(dp2lm),
         fm=expectdp2s,rfm=dmo[1],qfm=dmo[2],sdfm=sd(dp2m)
  )
  return(all)
}
rrayleigh<-function (n, scale = 1) {
  sample1 <- scale * sqrt(-2 * log(runif(n)))
  sample1[scale <= 0] <- NaN
  sample1
}

NRSs<-function (x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,standarddistribution=c("exponential","rayleigh"),SE=TRUE,SD=FALSE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  lengthx<-length(sortedx)
  if(standarddistribution=="exponential"){
    drm=0.3665
    dqm=0.82224
    
    dlrmscale=0.3659
    drmscale=0.7930
    dlqmscale=0.8218
    dqmscale=0.7825
    
    dlrmtm=0.1808
    drmtm=1.7492
    dlqmtm=1.1753
    dqmtm=0.5715
    
    dlrmfm=-0.3542
    drmfm=3.4560
    dlqmfm=NaN
    dqmfm=0.1246
    
  }else if (standarddistribution=="rayleigh"){
    drm=0.4025526
    dqm=0.4452798
    
    dlrmscale=0.3095063
    drmscale=0.3862421
    dlqmscale=0.7055909
    dqmscale=1.097661
    
    dlrmtm=0.1561753
    drmtm=0.7855876
    dlqmtm=0.8038741
    dqmtm=0.9621217
    
    dlrmfm=0.2290345
    drmfm=0.8758908
    dlqmfm=0.6571804
    dqmfm=0.7304692
  }
  mmm1<-mmme(x=sortedx,interval=9,fast=TRUE,batch=1000,sorted=TRUE,drm=drm,dqm=dqm)
  rqscale1<-rqscale(x=sortedx,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
  rqtm1<-rqtm(x=sortedx,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
  rqfm1<-rqfm(x=sortedx,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
  
  if(SD){
    all<-c(c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4]),
           c(l2=rqscale1[1],rl2=rqscale1[2],ql2=rqscale1[3],sdrql2=rqscale1[4]),
           c(sd=sqrt(rqscale1[5]),
             rsd=sqrt(rqscale1[6]),qsd=sqrt(rqscale1[7]),sdrsd=sqrt(rqscale1[6])*(1/2)*(rqscale1[8]/rqscale1[6]),sdqsd=sqrt(rqscale1[7])*(1/2)*(rqscale1[8]/rqscale1[7])
           ),
           c(l3=rqtm1[1]/rqscale1[1],rl3=rqtm1[2]/rqscale1[2],ql3=rqtm1[3]/rqscale1[3],sdrl3=(rqtm1[2]/rqscale1[2])*((rqtm1[4]/rqtm1[2])^2+(rqscale1[4]/rqscale1[2])^2)^(1/2),
             sdql3=(rqtm1[3]/rqscale1[3])*((rqtm1[4]/rqtm1[3])^2+(rqscale1[4]/rqscale1[3])^2)^(1/2)),
           
           c(skew=(rqtm1[5])/((rqscale1[5])^(3/2)),rskew=(rqtm1[6])/((rqscale1[6])^(3/2)),qskew=(rqtm1[7])/((rqscale1[7])^(3/2)),
             sdrskew=(rqtm1[6])/((rqscale1[6])^(3/2))*((rqtm1[8]/rqtm1[6])^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2)^(1/2),
             sdqskew=(rqtm1[7])/((rqscale1[7])^(3/2))*((rqtm1[8]/rqtm1[7])^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2)^(1/2)
           ),
           
           c(l4=rqfm1[1]/rqscale1[1],rl4=rqfm1[2]/rqscale1[2],ql4=rqfm1[3]/rqscale1[3],sdrl4=(rqfm1[2]/rqscale1[2])*((rqfm1[4]/rqfm1[2])^2+(rqscale1[4]/rqscale1[2])^2)^(1/2),
             sdql4=(rqfm1[3]/rqscale1[3])*((rqfm1[4]/rqfm1[3])^2+(rqscale1[4]/rqscale1[3])^2)^(1/2)
           ),
           
           c(kurt=(rqfm1[5])/((rqscale1[5])^(2)),rkurt=(rqfm1[6])/((rqscale1[6])^(2)),qkurt=(rqfm1[7])/((rqscale1[7])^(2)),
             sdrkurt=((rqfm1[6])/((rqscale1[6])^(2)))*((rqfm1[8]/rqfm1[6])^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2)^(1/2),
             sdqkurt=((rqfm1[7])/((rqscale1[7])^(2)))*((rqfm1[8]/rqfm1[7])^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2)^(1/2)
           ))
    
    
  }else if(SE){
    
    all<-c(c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4]),semean=sqrt(rqscale1[5])/sqrt(lengthx),
           c(l2=rqscale1[1],rl2=rqscale1[2],ql2=rqscale1[3],serql2=rqscale1[4]/sqrt(lengthx)),
           c(sd=sqrt(rqscale1[5]),
             rsd=sqrt(rqscale1[6]),qsd=sqrt(rqscale1[7]),sersd=sqrt(rqscale1[6])*(1/2)*(rqscale1[8]/rqscale1[6])/sqrt(lengthx),seqsd=sqrt(rqscale1[7])*(1/2)*(rqscale1[8]/rqscale1[7])/sqrt(lengthx)
           ),
           c(l3=rqtm1[1]/rqscale1[1],rl3=rqtm1[2]/rqscale1[2],ql3=rqtm1[3]/rqscale1[3],serl3=(1/sqrt(lengthx))*(rqtm1[2]/rqscale1[2])*((rqtm1[4]/rqtm1[2])^2+(rqscale1[4]/rqscale1[2])^2)^(1/2),
             seql3=(1/sqrt(lengthx))*(rqtm1[3]/rqscale1[3])*((rqtm1[4]/rqtm1[3])^2+(rqscale1[4]/rqscale1[3])^2)^(1/2)),
           
           c(skew=(rqtm1[5])/((rqscale1[5])^(3/2)),rskew=(rqtm1[6])/((rqscale1[6])^(3/2)),qskew=(rqtm1[7])/((rqscale1[7])^(3/2)),
             serskew=(1/sqrt(lengthx))*(rqtm1[6])/((rqscale1[6])^(3/2))*((rqtm1[8]/rqtm1[6])^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2)^(1/2),
             seqskew=(1/sqrt(lengthx))*(rqtm1[7])/((rqscale1[7])^(3/2))*((rqtm1[8]/rqtm1[7])^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2)^(1/2)
           ),
           
           c(l4=rqfm1[1]/rqscale1[1],rl4=rqfm1[2]/rqscale1[2],ql4=rqfm1[3]/rqscale1[3],serl4=(1/sqrt(lengthx))*(rqfm1[2]/rqscale1[2])*((rqfm1[4]/rqfm1[2])^2+(rqscale1[4]/rqscale1[2])^2)^(1/2),
             seql4=(1/sqrt(lengthx))*(rqfm1[3]/rqscale1[3])*((rqfm1[4]/rqfm1[3])^2+(rqscale1[4]/rqscale1[3])^2)^(1/2)
           ),
           
           c(kurt=(rqfm1[5])/((rqscale1[5])^(2)),rkurt=(rqfm1[6])/((rqscale1[6])^(2)),qkurt=(rqfm1[7])/((rqscale1[7])^(2)),
             serkurt=(1/sqrt(lengthx))*((rqfm1[6])/((rqscale1[6])^(2)))*((rqfm1[8]/rqfm1[6])^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2)^(1/2),
             seqkurt=(1/sqrt(lengthx))*((rqfm1[7])/((rqscale1[7])^(2)))*((rqfm1[8]/rqfm1[7])^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2)^(1/2)
           ))
  }
  
  if((rqfm1[7])/((rqscale1[7])^(2))<5.5 & standarddistribution=="exponential"){
    print("the kurtosis is too low, it might be better to use the rayleigh distribution as the standard distribution")
  }
  
  return(all)
}



#test
x<-rexp(5400,1)

#the population mean is 1
#the population variance is 1
#the population L2-moment is 1/2
#the population skewness is 2
#the population L-skewness is 1/3
#the population kurtosis is 9
#the population L-kurtosis is 1/6
#no d for ql4, because the distribution of U-statistic of L4-moment does not follow mean-ETM-median inequality

NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="exponential",SE=TRUE,SD=FALSE)

#sample mean, standard deviation, skewness, kurtosis are provided to compared,
#the standard deviations of the underlying distribution are used to estimate the standard error
#, for example, sdrkurt is ~90, a rough estimation of SE is 90/sqrt(5400)=1.224745
#if want to see the standard deviation of the distribution of U-statistic, instead of SE, just set SD to TRUE
NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="exponential",SE=FALSE,SD=TRUE)


x<-rgamma(5400,4,1)
#the population mean is 4
#the population variance is 4
#the population L2-moment is 1.09375
#the population skewness is 1
#the population L-skewness is 0.1646599
#the population kurtosis is 4.5
#the population L-kurtosis is 0.1312521
NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="exponential",SE=TRUE,SD=FALSE)
#a rule of thumb is for desired consistency performance (all four moments >90%) is the kurtosis of the underlying distribution 
#should be within [1/2,2] times that of the standard distribution used to calibrated the d values
#that means, using exponential as the standard distribution, the kurtosis should be within 4.5 to 18

#while accurately estimating population kurtosis is hard, find a rough range and choose the right standard should be easy in practice, 
#also, if the quantile kurtosis is less than 5, highly indicating the need to change to rayleigh
NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="rayleigh",SE=TRUE,SD=FALSE)

x<-rrayleigh(n=5400, scale = 1) 
NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="exponential",SE=TRUE,SD=FALSE)
NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="rayleigh",SE=TRUE,SD=FALSE)


x<-rpois(5400,8)
#the population mean is 8
#quantile mean returns 7 or 9
#the population variance is 8
#the population L2-moment is 1.583
#the population skewness is 0.3535534
#the population L-skewness is 0.0592
#the population kurtosis is 1/8+3
#the population L-kurtosis is 0.1204

#quantile L-moments is not suitable for discrete unimodal distributions

#because the kurtosis of poisson is close to 3, the biases of robust/quantile kurtosis based on the exponential distribution are large.

NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="exponential",SE=TRUE,SD=FALSE)

NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="rayleigh",SE=TRUE,SD=FALSE)


x<-c(rnorm(5400))

NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="exponential",SE=TRUE,SD=FALSE)
NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="rayleigh",SE=TRUE,SD=FALSE)

x<-c(rlogis(5400, location = 0, scale = 1))

NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="exponential",SE=TRUE,SD=FALSE)
NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="rayleigh",SE=TRUE,SD=FALSE)

#NRSs have excellent performance even for heavy tailed distributions. 
#Even the kurtosis is extreme or infinite, while the consistency is poor, it is still better than all current robust statistics.
a=500
library(VGAM)
x<-c(VGAM::rpareto(5400, scale  = 1, shape=2+a/100))

NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="exponential",SE=TRUE,SD=FALSE)
NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="rayleigh",SE=TRUE,SD=FALSE)

