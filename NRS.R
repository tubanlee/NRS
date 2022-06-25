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
mmm<-function(expectboot,expecttrue,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE,drm=0.366439,dqm=0.8241268){
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
  listd<-c(expecttrue,expectboot,rm1,qm1)
  return(listd)
}

mmme<-function(x,interval=9,fast=TRUE,batch=1000,sorted=FALSE,drm=0.366439,dqm=0.8241268){
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
  listd<-c(etm1=etm1[2],rm1=rm1,qm1=qm1)
  return(listd)
}

rqscale<-function (x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE){
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
  dlmo<-mmm(expectboot=mean(dp2lm),expecttrue=expectdps,x=dp2lm,interval=9,fast=TRUE,batch=10000,sorted=FALSE,drm=0.3658755,dqm=0.8225767)
  dmo<-mmm(expectboot=mean(dp2m),expecttrue=expectdp2s,x=dp2m,interval=9,fast=TRUE,batch=10000,sorted=FALSE,drm=0.7929569,dqm=0.7843745)
  all<-c(rl2=dlmo[3],ql2=dlmo[4],sdl2=sd(dp2lm),
         rvar=dmo[3],qvar=dmo[4],sdvar=sd(dp2m)
  )
  return(all)
}

rqtm<-function (x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE){
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
  dlmo<-mmm(expectboot=mean(dp2lm),expecttrue=expectdps,x=dp2lm,interval=9,fast=fast,batch=batch,sorted=FALSE,drm=0.1810246,dqm=1.175519)
  dmo<-mmm(expectboot=mean(dp2m),expecttrue=expectdp2s,x=dp2m,interval=9,fast=fast,batch=batch,sorted=FALSE,drm=1.74967,dqm=0.5715816)
  
  rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  
  alltm<-c(rlskew=dlmo[3]/rqscale1[1],qlskew=dlmo[4]/rqscale1[2],sdl3=sd(dp2lm),
         rskew=dmo[3]/((rqscale1[4])^(3/2)),qskew=dmo[4]/((rqscale1[5])^(3/2)),sdtm=sd(dp2m)
  )
  
  all<-c(rqscale1,alltm)
  
  return(all)
}

rqfm<-function (x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE){
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
  dlmo<-mmm(expectboot=mean(dp2lm),expecttrue=expectdps,x=dp2lm,interval=9,fast=fast,batch=batch,sorted=FALSE,drm=-0.3567702,dqm=NaN)
  dmo<-mmm(expectboot=mean(dp2m),expecttrue=expectdp2s,x=dp2m,interval=9,fast=fast,batch=batch,sorted=FALSE,drm=3.454845,dqm=0.1249966)
  
  rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  
  allfm<-c(rlkurt=dlmo[3]/rqscale1[1],qlkurt=dlmo[4]/rqscale1[2],sdl4=sd(dp2lm),
           rkurt=dmo[3]/((rqscale1[4])^(2)),qkurt=dmo[4]/((rqscale1[5])^(2)),sdfm=sd(dp2m)
  )
  
  all<-c(allfm)
  return(all)
}



#test
x<-rexp(5400,1)
#the population mean is 1
mean(x)
mmme(x,interval=9,fast=TRUE,batch=1000,sorted=FALSE,drm=0.366439,dqm=0.8241268)
#the population variance is 1
#the population L2-moment is 1/2
#the population skewness is 2
#the population L-skewness is 1/3
rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
#the population kurtosis is 9
#the population L-kurtosis is 1/6
rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)

#no d for ql4, because the distribution of U-statistic of L4-moment does not follow mean-ETM-median inequality

x<-rgamma(5400,4,1)
#the population mean is 4
mean(x)
mmme(x,interval=9,fast=TRUE,batch=1000,sorted=FALSE,drm=0.366439,dqm=0.8241268)
#the population variance is 4
#the population L2-moment is 1.09375
#the population skewness is 1
#the population L-skewness is 0.1646599
rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
#the population kurtosis is 4.5
#the population L-kurtosis is 0.1312521
rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)



x<-rpois(5400,8)
#the population mean is 8
mean(x)
#quantile mean returns 7 or 9
mmme(x,interval=9,fast=TRUE,batch=1000,sorted=FALSE,drm=0.366439,dqm=0.8241268)
#the population variance is 8
#the population L2-moment is 1.583
#the population skewness is 0.3535534
#the population L-skewness is 0.0592
rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
#the population kurtosis is 1/8+3
#the population L-kurtosis is 0.1204
rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)

#quantile L-moments is not suitable for discrete unimodal distributions

#because the kurtosis of poisson is close to 3, the biases of robust/quantile kurtosis are large.

#for more tests, use the codes in consistency.R

#bootstrap is asymptotic valid over all the distributions tested.
