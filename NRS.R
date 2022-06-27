


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
rqmean<-function (x, interval=9,fast=TRUE,batch=1000,sorted=FALSE,drm=0.3665,dqm=0.82224,cise = FALSE,alpha = 0.05,nboot = 1000){
  if(cise){
    return (mmmeci(x, interval=interval,fast=fast,batch=batch,sorted=sorted,drm=drm,dqm=dqm,alpha = alpha,nboot = nboot))
  } else {return (mmme(x,interval=interval,fast=fast,batch=batch,sorted=sorted,drm=drm,dqm=dqm))}
}

mmmeci<-function(x,interval=9,fast=TRUE,batch=1000,sorted=FALSE,drm=0.3665,dqm=0.82224,alpha=0.05,nboot=1000){
  data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  pairs1 <- data.frame()
  pairs2 <- data.frame()
  pairs3 <- data.frame()
  pairs4 <- data.frame()
  for (i in 1:nboot) {
    robustlocation1<-mmme(data[i,],interval=interval,fast=fast,batch=batch,sorted=sorted,drm=drm,dqm=dqm)
    pairs1 <- rbind(pairs1,robustlocation1[1])
    pairs2 <- rbind(pairs2,robustlocation1[2])
    pairs3 <- rbind(pairs3,robustlocation1[3])
    pairs4 <- rbind(pairs4,robustlocation1[4])
  }
  bootlist1<-sort(as.matrix(pairs1))
  bootlist2<-sort(as.matrix(pairs2))
  bootlist3<-sort(as.matrix(pairs3))
  bootlist4<-sort(as.matrix(pairs4))
  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  estimate=mmme(x,interval=interval,fast=fast,batch=batch,sorted=sorted,drm=drm,dqm=dqm)
  result <- list(cimean=c(bootlist1[low],bootlist1[up]), cietm=c(bootlist2[low],bootlist2[up]),cirm=c(bootlist3[low],bootlist3[up]),
                 ciqm=c(bootlist4[low],bootlist4[up]),semean=sd(bootlist1),seetm=sd(bootlist2),serm=sd(bootlist3),seqm=sd(bootlist4),estimate=estimate)
  return(result)
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
rqscaleci<-function(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,dlrm=0.3659,drm=0.7930,dlqm=0.8218,dqm=0.7825,alpha=0.05,nboot=1000){
  data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  pairs1 <- data.frame()
  pairs2 <- data.frame()
  pairs3 <- data.frame()
  pairs4 <- data.frame()
  pairs5 <- data.frame()
  pairs6 <- data.frame()
  for (i in 1:nboot) {
    robustlocation1<-rqscale(data[i,],interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=sorted,dlrm=dlrm,drm=drm,dlqm=dlqm,dqm=dqm)
    pairs1 <- rbind(pairs1,robustlocation1[1])
    pairs2 <- rbind(pairs2,robustlocation1[2])
    pairs3 <- rbind(pairs3,robustlocation1[3])
    pairs4 <- rbind(pairs4,sqrt(robustlocation1[5]))
    pairs5 <- rbind(pairs5,sqrt(robustlocation1[6]))
    pairs6 <- rbind(pairs6,sqrt(robustlocation1[7]))
  }
  bootlist1<-sort(as.matrix(pairs1))
  bootlist2<-sort(as.matrix(pairs2))
  bootlist3<-sort(as.matrix(pairs3))
  bootlist4<-sort(as.matrix(pairs4))
  bootlist5<-sort(as.matrix(pairs5))
  bootlist6<-sort(as.matrix(pairs6))
  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  estimate=rqscale(x,interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=sorted,dlrm=dlrm,drm=drm,dlqm=dlqm,dqm=dqm)
  result <- list(cil2=c(bootlist1[low],bootlist1[up]), cirl2=c(bootlist2[low],bootlist2[up]),ciql2=c(bootlist3[low],bootlist3[up]),
                 civar=c(bootlist4[low],bootlist4[up]),cirvar=c(bootlist5[low],bootlist5[up]),ciqvar=c(bootlist6[low],bootlist6[up]),
                 sel2=sd(bootlist1),serl2=sd(bootlist2),seql2=sd(bootlist3),sevar=sd(bootlist4),servar=sd(bootlist5),seqvar=sd(bootlist6),
                 estimate=estimate
                 )
  return(result)
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
rqtmci<-function(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,dlrm=0.1808,drm=1.7492,dlqm=1.1753,dqm=0.5715,alpha=0.05,nboot=1000){
  data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  pairs1 <- data.frame()
  pairs2 <- data.frame()
  pairs3 <- data.frame()
  pairs4 <- data.frame()
  pairs5 <- data.frame()
  pairs6 <- data.frame()
  for (i in 1:nboot) {
    robustlocation1<-rqtm(data[i,],interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=sorted,dlrm=dlrm,drm=drm,dlqm=dlqm,dqm=dqm)
    pairs1 <- rbind(pairs1,robustlocation1[1])
    pairs2 <- rbind(pairs2,robustlocation1[2])
    pairs3 <- rbind(pairs3,robustlocation1[3])
    pairs4 <- rbind(pairs4,robustlocation1[5])
    pairs5 <- rbind(pairs5,robustlocation1[6])
    pairs6 <- rbind(pairs6,robustlocation1[7])
  }
  bootlist1<-sort(as.matrix(pairs1))
  bootlist2<-sort(as.matrix(pairs2))
  bootlist3<-sort(as.matrix(pairs3))
  bootlist4<-sort(as.matrix(pairs4))
  bootlist5<-sort(as.matrix(pairs5))
  bootlist6<-sort(as.matrix(pairs6))
  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  estimate=rqtm(x,interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=sorted,dlrm=dlrm,drm=drm,dlqm=dlqm,dqm=dqm)
  result <- list(cil3=c(bootlist1[low],bootlist1[up]), cirl3=c(bootlist2[low],bootlist2[up]),ciql3=c(bootlist3[low],bootlist3[up]),
                 citm=c(bootlist4[low],bootlist4[up]),cirtm=c(bootlist5[low],bootlist5[up]),ciqtm=c(bootlist6[low],bootlist6[up]),
                 sel3=sd(bootlist1),serl3=sd(bootlist2),seql3=sd(bootlist3),setm=sd(bootlist4),sertm=sd(bootlist5),seqtm=sd(bootlist6),
                 estimate=estimate
  )
  return(result)
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

rqfmci<-function(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,dlrm=-0.3542,drm=3.4560,dlqm=NaN,dqm=0.1246,alpha=0.05,nboot=1000){
  data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  pairs1 <- data.frame()
  pairs2 <- data.frame()
  pairs3 <- data.frame()
  pairs4 <- data.frame()
  pairs5 <- data.frame()
  pairs6 <- data.frame()
  for (i in 1:nboot) {
    robustlocation1<-rqfm(data[i,],interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=sorted,dlrm=dlrm,drm=drm,dlqm=dlqm,dqm=dqm)
    pairs1 <- rbind(pairs1,robustlocation1[1])
    pairs2 <- rbind(pairs2,robustlocation1[2])
    pairs3 <- rbind(pairs3,robustlocation1[3])
    pairs4 <- rbind(pairs4,robustlocation1[5])
    pairs5 <- rbind(pairs5,robustlocation1[6])
    pairs6 <- rbind(pairs6,robustlocation1[7])
  }
  bootlist1<-sort(as.matrix(pairs1))
  bootlist2<-sort(as.matrix(pairs2))
  bootlist3<-sort(as.matrix(pairs3))
  bootlist4<-sort(as.matrix(pairs4))
  bootlist5<-sort(as.matrix(pairs5))
  bootlist6<-sort(as.matrix(pairs6))
  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  estimate=rqfm(x,interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=sorted,dlrm=dlrm,drm=drm,dlqm=dlqm,dqm=dqm)
  result <- list(cil4=c(bootlist1[low],bootlist1[up]), cirl4=c(bootlist2[low],bootlist2[up]),ciql4=c(bootlist3[low],bootlist3[up]),
                 cifm=c(bootlist4[low],bootlist4[up]),cirfm=c(bootlist5[low],bootlist5[up]),ciqfm=c(bootlist6[low],bootlist6[up]),
                 sel4=sd(bootlist1),serl4=sd(bootlist2),seql4=sd(bootlist3),sefm=sd(bootlist4),serfm=sd(bootlist5),seqfm=sd(bootlist6),
                 estimate=estimate
  )
  return(result)
}

rrayleigh<-function (n, scale = 1) {
  sample1 <- scale * sqrt(-2 * log(runif(n)))
  sample1[scale <= 0] <- NaN
  sample1
}

NRSssimple<-function (x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,standarddistribution=c("exponential","rayleigh"),SE=TRUE,SD=FALSE){
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
  mmm1<-mmme(x=sortedx,interval=interval,fast=fast,batch=batch,sorted=sorted,drm=drm,dqm=dqm)
  rqscale1<-rqscale(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=TRUE,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
  rqtm1<-rqtm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=TRUE,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
  rqfm1<-rqfm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=TRUE,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
  
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
             rsd=sqrt(rqscale1[6]),qsd=sqrt(rqscale1[7]),sersd=sqrt(rqscale1[6])*(1/2)*((rqscale1[8])/rqscale1[6])/sqrt(lengthx),seqsd=sqrt(rqscale1[7])*(1/2)*(rqscale1[8]/rqscale1[7])/sqrt(lengthx)
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
NRSs<-function(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,standarddistribution=c("exponential","rayleigh"),SE=TRUE,SD=FALSE,cise = FALSE,alpha = 0.05,nboot = 100){
  if(cise){
    return (NRSsci(x, interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=sorted,standarddistribution=standarddistribution,alpha=alpha,nboot=nboot))
  } else {return (NRSssimple(x, interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=sorted,standarddistribution=standarddistribution,SE=SE,SD = SD))
    
}}
NRSsci<-function (x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,standarddistribution=c("exponential","rayleigh"),alpha=0.05,nboot=100){
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
  data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  pairs1 <- data.frame()
  pairs2 <- data.frame()
  pairs3 <- data.frame()
  pairs4 <- data.frame()
  pairs5 <- data.frame()
  pairs6 <- data.frame()
  pairs7 <- data.frame()
  pairs8 <- data.frame()
  pairs9 <- data.frame()
  pairs10 <- data.frame()
  pairs11 <- data.frame()
  pairs12 <- data.frame()
  pairs13 <- data.frame()
  pairs14 <- data.frame()
  pairs15 <- data.frame()
  pairs16 <- data.frame()
  pairs17 <- data.frame()
  pairs18 <- data.frame()
  pairs19 <- data.frame()
  pairs20 <- data.frame()
  pairs21 <- data.frame()
  pairs22 <- data.frame()
  for (i in 1:nboot) {
    mmm1<-mmme(x=data[i,],interval=interval,fast=fast,batch=batch,sorted=sorted,drm=drm,dqm=dqm)
    rqscale1<-rqscale(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=TRUE,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
    rqtm1<-rqtm(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=TRUE,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
    rqfm1<-rqfm(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=TRUE,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
    
    estimate1<-c(c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4]),
                c(l2=rqscale1[1],rl2=rqscale1[2],ql2=rqscale1[3]),
                c(sd=sqrt(rqscale1[5]),
                  rsd=sqrt(rqscale1[6]),qsd=sqrt(rqscale1[7])
                ),
                c(l3=rqtm1[1]/rqscale1[1],rl3=rqtm1[2]/rqscale1[2],ql3=rqtm1[3]/rqscale1[3]),
                
                c(skew=(rqtm1[5])/((rqscale1[5])^(3/2)),rskew=(rqtm1[6])/((rqscale1[6])^(3/2)),qskew=(rqtm1[7])/((rqscale1[7])^(3/2))),
                
                c(l4=rqfm1[1]/rqscale1[1],rl4=rqfm1[2]/rqscale1[2],ql4=rqfm1[3]/rqscale1[3]
                ),
                
                c(kurt=(rqfm1[5])/((rqscale1[5])^(2)),rkurt=(rqfm1[6])/((rqscale1[6])^(2)),qkurt=(rqfm1[7])/((rqscale1[7])^(2))))
    pairs1 <- rbind(pairs1,estimate1[1])
    pairs2 <- rbind(pairs2,estimate1[2])
    pairs3 <- rbind(pairs3,estimate1[3])
    pairs4 <- rbind(pairs4,estimate1[4])
    pairs5 <- rbind(pairs5,estimate1[5])
    pairs6 <- rbind(pairs6,estimate1[6])
    pairs7 <- rbind(pairs7,estimate1[7])
    pairs8 <- rbind(pairs8,estimate1[8])
    pairs9 <- rbind(pairs9,estimate1[9])
    pairs10 <- rbind(pairs10,estimate1[10])
    pairs11 <- rbind(pairs11,estimate1[11])
    pairs12 <- rbind(pairs12,estimate1[12])
    pairs13 <- rbind(pairs13,estimate1[13])
    pairs14 <- rbind(pairs14,estimate1[14])
    pairs15 <- rbind(pairs15,estimate1[15])
    pairs16 <- rbind(pairs16,estimate1[16])
    pairs17 <- rbind(pairs17,estimate1[17])
    pairs18 <- rbind(pairs18,estimate1[18])
    pairs19 <- rbind(pairs19,estimate1[19])
    pairs20 <- rbind(pairs20,estimate1[20])
    pairs21 <- rbind(pairs21,estimate1[21])
    pairs22 <- rbind(pairs22,estimate1[22])
  }
  bootlist1<-sort(as.matrix(pairs1))
  bootlist2<-sort(as.matrix(pairs2))
  bootlist3<-sort(as.matrix(pairs3))
  bootlist4<-sort(as.matrix(pairs4))
  bootlist5<-sort(as.matrix(pairs5))
  bootlist6<-sort(as.matrix(pairs6))
  bootlist7<-sort(as.matrix(pairs7))
  bootlist8<-sort(as.matrix(pairs8))
  bootlist9<-sort(as.matrix(pairs9))
  bootlist10<-sort(as.matrix(pairs10))
  bootlist11<-sort(as.matrix(pairs11))
  bootlist12<-sort(as.matrix(pairs12))
  bootlist13<-sort(as.matrix(pairs13))
  bootlist14<-sort(as.matrix(pairs14))
  bootlist15<-sort(as.matrix(pairs15))
  bootlist16<-sort(as.matrix(pairs16))
  bootlist17<-sort(as.matrix(pairs17))
  bootlist18<-sort(as.matrix(pairs18))
  bootlist19<-sort(as.matrix(pairs19))
  bootlist20<-sort(as.matrix(pairs20))
  bootlist21<-sort(as.matrix(pairs21))
  bootlist22<-sort(as.matrix(pairs22))
  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  mmm1<-mmme(x=sortedx,interval=interval,fast=fast,batch=batch,sorted=sorted,drm=drm,dqm=dqm)
  rqscale1<-rqscale(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=TRUE,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
  rqtm1<-rqtm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=TRUE,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
  rqfm1<-rqfm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,subsample=subsample,sorted=TRUE,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
  estimate<-c(c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4]),
           c(l2=rqscale1[1],rl2=rqscale1[2],ql2=rqscale1[3]),
           c(sd=sqrt(rqscale1[5]),
             rsd=sqrt(rqscale1[6]),qsd=sqrt(rqscale1[7])
           ),
           c(l3=rqtm1[1]/rqscale1[1],rl3=rqtm1[2]/rqscale1[2],ql3=rqtm1[3]/rqscale1[3]),
           
           c(skew=(rqtm1[5])/((rqscale1[5])^(3/2)),rskew=(rqtm1[6])/((rqscale1[6])^(3/2)),qskew=(rqtm1[7])/((rqscale1[7])^(3/2))),
           
           c(l4=rqfm1[1]/rqscale1[1],rl4=rqfm1[2]/rqscale1[2],ql4=rqfm1[3]/rqscale1[3]
           ),
           
           c(kurt=(rqfm1[5])/((rqscale1[5])^(2)),rkurt=(rqfm1[6])/((rqscale1[6])^(2)),qkurt=(rqfm1[7])/((rqscale1[7])^(2))))
  ci<-c(c(mean=c(bootlist1[low],bootlist1[up]),etm=c(bootlist2[low],bootlist2[up]),rm=c(bootlist3[low],bootlist3[up]),qm=c(bootlist4[low],bootlist4[up])),
              c(l2=c(bootlist5[low],bootlist5[up]),rl2=c(bootlist6[low],bootlist6[up]),ql2=c(bootlist7[low],bootlist7[up])),
              c(sd=c(bootlist8[low],bootlist8[up]),
                rsd=c(bootlist9[low],bootlist9[up]),qsd=c(bootlist10[low],bootlist10[up])
              ),
              c(l3=c(bootlist11[low],bootlist11[up]),rl3=c(bootlist12[low],bootlist12[up]),ql3=c(bootlist13[low],bootlist13[up])),
              
              c(skew=c(bootlist14[low],bootlist14[up]),rskew=c(bootlist15[low],bootlist15[up]),qskew=c(bootlist16[low],bootlist16[up])),
              
              c(l4=c(bootlist17[low],bootlist17[up]),rl4=c(bootlist18[low],bootlist18[up]),ql4=c(bootlist19[low],bootlist19[up])
              ),
              
              c(kurt=c(bootlist20[low],bootlist20[up]),rkurt=c(bootlist21[low],bootlist21[up]),qkurt=c(bootlist22[low],bootlist2[up])))
  
  se<-c(c(mean=sd(bootlist1),etm=sd(bootlist2),rm=sd(bootlist3),qm=sd(bootlist4)),
        c(l2=sd(bootlist5),rl2=sd(bootlist6),ql2=sd(bootlist7)),
        c(sd=sd(bootlist8),
          rsd=sd(bootlist9),qsd=sd(bootlist10)
        ),
        c(l3=sd(bootlist11),rl3=sd(bootlist12),ql3=sd(bootlist13)),
        
        c(skew=sd(bootlist14),rskew=sd(bootlist15),qskew=sd(bootlist16)),
        
        c(l4=sd(bootlist17),rl4=sd(bootlist18),ql4=sd(bootlist19)
        ),
        
        c(kurt=sd(bootlist20),rkurt=sd(bootlist21),qkurt=sd(bootlist22)))
  
  all<-c(ci=ci,se=se,estimate=estimate)
  
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
#if want to see the standard deviation of the distribution of U-statistic,  set SD to TRUE.
#this standard deviation is calculated based on the law of propogatio of uncertainty.
NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,standarddistribution="exponential",SE=FALSE,SD=TRUE,cise = FALSE,alpha = 0.05,nboot = 100)

#sample mean, standard deviation, skewness, kurtosis are provided to compared,
#the standard deviations of the underlying distribution can be used to estimate the standard error
#, for example, sdrkurt is ~90, a rough estimation of SE is 90/sqrt(5400)=1.224745 
#but this is just rough estimation, since the influence functions haven't been derived yet. 

NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,standarddistribution="exponential",SE=TRUE,SD=FALSE,cise = FALSE,alpha = 0.05,nboot = 100)
#the standard error and confidential interval of robust/quantile mean can also be estimated based on bootstrap

rqmean(x, interval=9,fast=TRUE,batch=1000,sorted=FALSE,drm=0.3665,dqm=0.82224,cise = TRUE,alpha = 0.05,nboot = 1000)
#similar approach can be applied to all NRSs, but the challenge is the computational time (just 100 nboot takes ~10 mins)

NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE,standarddistribution="exponential",cise = TRUE,alpha = 0.05,nboot = 100)

#comparing the above results implies that the rough estimation is generally well (the average deviation is ~30%)


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

#quantile L-moments are not suitable for discrete unimodal distributions

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

a=500
library(VGAM)
x<-c(VGAM::rpareto(5400, scale  = 1, shape=2+a/100))

NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="exponential",SE=TRUE,SD=FALSE)
NRSs(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=FALSE,standarddistribution="rayleigh",SE=TRUE,SD=FALSE)



#for more tests, use the codes in consistency.R
