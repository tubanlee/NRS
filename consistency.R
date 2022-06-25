



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
  names(qm1)<-NULL
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
  names(qm1)<-NULL
  rm1<--drm*etm1[3]+etm1[2]+drm*etm1[2]
  listd<-c(etm1[2],rm1,qm1)
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
  all<-c(rl3=dlmo[3],ql3=dlmo[4],sdl3=sd(dp2lm),
         rtm=dmo[3],qtm=dmo[4],sdtm=sd(dp2m)
  )
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
  all<-c(rl4=dlmo[3],ql4=dlmo[4],sdl4=sd(dp2lm),
         rfm=dmo[3],qfm=dmo[4],sdfm=sd(dp2m)
  )
  return(all)
}

winsor<-function (x, fraction=.05)
{
  lim <- quantile(x, probs=c(fraction, 1-fraction))
  x[ x < lim[1] ] <- lim[1]
  x[ x > lim[2] ] <- lim[2]
  x
}
library(Lmoments)

#the biases of NRSs for laplace are very small, because the kurtosis is 6, not much differ from that of exponential distribution (9)
allforlaplace<-c()
for(i in (1:10)){
  library(VGAM)
  x<-c(rlaplace(5400, location = 0, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  mean1<-mean(x)
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
  rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  #rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=FALSE,subsample=54000,sorted=TRUE)
  rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  targetm<-0
  targetvar<-(2)
  targettm<-0
  targetfm<-6*(sqrt(2)^4)
  targetl2<-3/4
  targetl3<-0
  targetl4<-(3/4)*1/(3*sqrt(2))
  all<-c(
    meanBias=((mean1-targetm)/(sd1)),
    tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1[3]-targetm)/(sd1)),rmBias=((rqmmm1[2]-targetm)/(sd1)),qmBias=((rqmmm1[3]-targetm)/(sd1)),
    rvarBias=((rqscale1[4]-targetvar)/(rqscale1[6])),qvarBias=((rqscale1[5]-targetvar)/(rqscale1[6])),rl2Bias=((rqscale1[1]-targetl2)/(rqscale1[3])),ql2Bias=((rqscale1[2]-targetl2)/(rqscale1[3])),
    rtmBias=((rqtm1[4]-targettm)/(rqtm1[6])),qtmBias=((rqtm1[5]-targettm)/(rqtm1[6])),rl3Bias=((rqtm1[1]-targetl3)/(rqtm1[3])),ql3Bias=((rqtm1[2]-targetl3)/(rqtm1[3])),
    rfmBias=(rqfm1[4]-targetfm)/(rqfm1[6]),qfmBias=(rqfm1[5]-targetfm)/(rqfm1[6]),rl4Bias=((rqfm1[1]-targetl4)/(rqfm1[3])),ql4Bias=((rqfm1[2]-targetl4)/(rqfm1[3])))
  allforlaplace<-rbind(allforlaplace,all)
}
allforlaplace[is.infinite(allforlaplace)] <-NA

write.csv(allforlaplace,paste("Conlaplace,",batchnumber,".csv", sep = ","), row.names = TRUE)


#the biases of robust/quantile fourth moments for normal are huge, because the kurtosis is 3, too differ from that of exponential distribution (9)

#but robust 4th moment still remain high consistent, because the distributions of U-statistics of L-moments are more symmetric.

allfornorm<-c()
for(i in (1:10)){
  x<-c(rnorm(5400))
  x<-sort(x,decreasing = FALSE,method ="radix")
  mean1<-mean(x)
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-mmme(x,interval=9,fast=TRUE,batch=1000,sorted=FALSE,drm=0.366439,dqm=0.8241268)
  rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  #rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=FALSE,subsample=54000,sorted=TRUE)
  rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  targetm<-0
  targetvar<-1
  targettm<-0
  targetfm<-3
  targetl2<-1/sqrt(pi)
  targetl3<-0
  targetl4<-(1/sqrt(pi))*(30*(1/(pi))*(atan(sqrt(2)))-9)
  all<-c(
    meanBias=((mean1-targetm)/(sd1)),
    tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1[3]-targetm)/(sd1)),rmBias=((rqmmm1[2]-targetm)/(sd1)),qmBias=((rqmmm1[3]-targetm)/(sd1)),
    rvarBias=((rqscale1[4]-targetvar)/(rqscale1[6])),qvarBias=((rqscale1[5]-targetvar)/(rqscale1[6])),rl2Bias=((rqscale1[1]-targetl2)/(rqscale1[3])),ql2Bias=((rqscale1[2]-targetl2)/(rqscale1[3])),
    rtmBias=((rqtm1[4]-targettm)/(rqtm1[6])),qtmBias=((rqtm1[5]-targettm)/(rqtm1[6])),rl3Bias=((rqtm1[1]-targetl3)/(rqtm1[3])),ql3Bias=((rqtm1[2]-targetl3)/(rqtm1[3])),
    rfmBias=(rqfm1[4]-targetfm)/(rqfm1[6]),qfmBias=(rqfm1[5]-targetfm)/(rqfm1[6]),rl4Bias=((rqfm1[1]-targetl4)/(rqfm1[3])),ql4Bias=((rqfm1[2]-targetl4)/(rqfm1[3])))
  allfornorm<-rbind(allfornorm,all)
}
allfornorm[is.infinite(allfornorm)] <-NA

write.csv(allfornorm,paste("Connorm,",batchnumber,".csv", sep = ","), row.names = TRUE)

#the biases of robust/quantile fourth moments for logis are not very large ~0.2 and ~0.06, because the kurtosis is 4.2, better than normal

allforlogis<-c()
for(i in (1:10)){
  x<-c(rlogis(5400, location = 0, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  mean1<-mean(x)
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
  rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  #rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=FALSE,subsample=54000,sorted=TRUE)
  rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  targetm<-0
  targetvar<-((pi^2)/3)
  targettm<-0
  targetfm<-((6/5)+3)*((sqrt((pi^2)/3))^4)
  targetl2<-1
  targetl3<-0
  targetl4<-1/6
  all<-c(
    meanBias=((mean1-targetm)/(sd1)),
    tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1[3]-targetm)/(sd1)),rmBias=((rqmmm1[2]-targetm)/(sd1)),qmBias=((rqmmm1[3]-targetm)/(sd1)),
    rvarBias=((rqscale1[4]-targetvar)/(rqscale1[6])),qvarBias=((rqscale1[5]-targetvar)/(rqscale1[6])),rl2Bias=((rqscale1[1]-targetl2)/(rqscale1[3])),ql2Bias=((rqscale1[2]-targetl2)/(rqscale1[3])),
    rtmBias=((rqtm1[4]-targettm)/(rqtm1[6])),qtmBias=((rqtm1[5]-targettm)/(rqtm1[6])),rl3Bias=((rqtm1[1]-targetl3)/(rqtm1[3])),ql3Bias=((rqtm1[2]-targetl3)/(rqtm1[3])),
    rfmBias=(rqfm1[4]-targetfm)/(rqfm1[6]),qfmBias=(rqfm1[5]-targetfm)/(rqfm1[6]),rl4Bias=((rqfm1[1]-targetl4)/(rqfm1[3])),ql4Bias=((rqfm1[2]-targetl4)/(rqfm1[3])))
  allforlogis<-rbind(allforlogis,all)
}
allforlogis[is.infinite(allforlogis)] <-NA

write.csv(allforlogis,paste("Conlogis,",batchnumber,".csv", sep = ","), row.names = TRUE)

#the biases of robust/quantile fourth moments for Rayleigh are very large ~0.3 and ~0.11, because the kurtosis is 3.245, very close to normal

allforRayleigh<-c()
for(i in (1:10)){
  library(VGAM)
  x<-c(rrayleigh(5400, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  mean1<-mean(x)
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
  rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  #rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=FALSE,subsample=54000,sorted=TRUE)
  rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  targetm<-sqrt(pi/2)
  targetvar<-(2-(pi/2))
  targettm<-((sqrt(2-(pi/2)))^3)*2*sqrt((pi))*(pi-3)/((4-pi)^(3/2))
  targetfm<-((sqrt(2-(pi/2)))^4)*(3-(6*(pi)^2-24*(pi)+16)/((4-pi)^(2)))
  targetl2<-0.5*(sqrt(2)-1)*sqrt(pi)
  targetl3<-(1/6)*(2*sqrt(6)+3*sqrt(2)-9)*sqrt(pi)
  targetl4<-(sqrt((77/6)-5*sqrt(6))-3/4)*sqrt(2*pi)
  all<-c(
    meanBias=((mean1-targetm)/(sd1)),
    tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1[3]-targetm)/(sd1)),rmBias=((rqmmm1[2]-targetm)/(sd1)),qmBias=((rqmmm1[3]-targetm)/(sd1)),
    rvarBias=((rqscale1[4]-targetvar)/(rqscale1[6])),qvarBias=((rqscale1[5]-targetvar)/(rqscale1[6])),rl2Bias=((rqscale1[1]-targetl2)/(rqscale1[3])),ql2Bias=((rqscale1[2]-targetl2)/(rqscale1[3])),
    rtmBias=((rqtm1[4]-targettm)/(rqtm1[6])),qtmBias=((rqtm1[5]-targettm)/(rqtm1[6])),rl3Bias=((rqtm1[1]-targetl3)/(rqtm1[3])),ql3Bias=((rqtm1[2]-targetl3)/(rqtm1[3])),
    rfmBias=(rqfm1[4]-targetfm)/(rqfm1[6]),qfmBias=(rqfm1[5]-targetfm)/(rqfm1[6]),rl4Bias=((rqfm1[1]-targetl4)/(rqfm1[3])),ql4Bias=((rqfm1[2]-targetl4)/(rqfm1[3])))
  allforRayleigh<-rbind(allforRayleigh,all)
}
allforRayleigh[is.infinite(allforRayleigh)] <-NA

write.csv(allforRayleigh,paste("ConRayleigh,",batchnumber,".csv", sep = ","), row.names = TRUE)


allforexp<-c()
for(i in (1:10)){
  x<-c(rexp(5400,1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  mean1<-mean(x)
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
  rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  #rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=FALSE,subsample=54000,sorted=TRUE)
  rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  targetm<-1
  targetvar<-1
  targettm<-2
  targetfm<-9
  targetl2<-1/2
  targetl3<-1/6
  targetl4<-1/12
  all<-c(
    meanBias=((mean1-targetm)/(sd1)),
    tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1[3]-targetm)/(sd1)),rmBias=((rqmmm1[2]-targetm)/(sd1)),qmBias=((rqmmm1[3]-targetm)/(sd1)),
    rvarBias=((rqscale1[4]-targetvar)/(rqscale1[6])),qvarBias=((rqscale1[5]-targetvar)/(rqscale1[6])),rl2Bias=((rqscale1[1]-targetl2)/(rqscale1[3])),ql2Bias=((rqscale1[2]-targetl2)/(rqscale1[3])),
    rtmBias=((rqtm1[4]-targettm)/(rqtm1[6])),qtmBias=((rqtm1[5]-targettm)/(rqtm1[6])),rl3Bias=((rqtm1[1]-targetl3)/(rqtm1[3])),ql3Bias=((rqtm1[2]-targetl3)/(rqtm1[3])),
    rfmBias=(rqfm1[4]-targetfm)/(rqfm1[6]),qfmBias=(rqfm1[5]-targetfm)/(rqfm1[6]),rl4Bias=((rqfm1[1]-targetl4)/(rqfm1[3])),ql4Bias=((rqfm1[2]-targetl4)/(rqfm1[3])))
  allforexp<-rbind(allforexp,all)
}
allforexp[is.infinite(allforexp)] <-NA

write.csv(allforexp,paste("Conexp,",batchnumber,".csv", sep = ","), row.names = TRUE)


#same issue, when kurtosis lower than 4, the performances of robust/quantile fourth moments are very poor.


#the 0.01-0.1 shape parameter range were removed, as the variance too high, often introducing bugs.
library(lmom)
listWeibull<-data.frame()
for (a in (10:300)) {
  allforweibull<-c()
  for(i in (1:1)){
    x<-c(rweibull(5400, shape=a/100, scale = 1))
    targetwei<-lmrwei(para = c(0, 1, a/100), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    mean1<-mean(x)
    tm1<-mean(x,1/9)
    wm1<-mean(winsor(x,fraction=1/9))
    sd1<-sd(x) 
    rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
    rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    #rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=FALSE,subsample=54000,sorted=TRUE)
    rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    targetm<-gamma(1+1/(a/100))
    targetvar<-(gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2)
    targettm<-((sqrt(gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2))^3)*(gamma(1+3/(a/100))-3*(gamma(1+1/(a/100)))*((gamma(1+2/(a/100))))+2*((gamma(1+1/(a/100)))^3))/((sqrt(gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2))^(3))
    targetfm<-((sqrt(gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2))^4)*(gamma(1+4/(a/100))-4*(gamma(1+3/(a/100)))*((gamma(1+1/(a/100))))+6*(gamma(1+2/(a/100)))*((gamma(1+1/(a/100)))^2)-3*((gamma(1+1/(a/100)))^4))/(((gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2))^(2))
    targetl2<-targetwei[2]
    targetl3<-targetwei[3]*targetwei[2]
    targetl4<-targetwei[4]*targetwei[2]
    all<-c(
      meanBias=((mean1-targetm)/(sd1)),
      tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1[3]-targetm)/(sd1)),rmBias=((rqmmm1[2]-targetm)/(sd1)),qmBias=((rqmmm1[3]-targetm)/(sd1)),
      rvarBias=((rqscale1[4]-targetvar)/(rqscale1[6])),qvarBias=((rqscale1[5]-targetvar)/(rqscale1[6])),rl2Bias=((rqscale1[1]-targetl2)/(rqscale1[3])),ql2Bias=((rqscale1[2]-targetl2)/(rqscale1[3])),
      rtmBias=((rqtm1[4]-targettm)/(rqtm1[6])),qtmBias=((rqtm1[5]-targettm)/(rqtm1[6])),rl3Bias=((rqtm1[1]-targetl3)/(rqtm1[3])),ql3Bias=((rqtm1[2]-targetl3)/(rqtm1[3])),
      rfmBias=(rqfm1[4]-targetfm)/(rqfm1[6]),qfmBias=(rqfm1[5]-targetfm)/(rqfm1[6]),rl4Bias=((rqfm1[1]-targetl4)/(rqfm1[3])),ql4Bias=((rqfm1[2]-targetl4)/(rqfm1[3])))
    
    allforweibull<-rbind(allforweibull,all)
  }
  allforweibull[is.infinite(allforweibull)] <-NA
  listWeibull<-rbind(listWeibull,allforweibull)
}

write.csv(listWeibull,paste("ConWeibull",batchnumber,".csv", sep = ","), row.names = TRUE)


library(lmom)
listgamma<-data.frame()
for (a in (1:300)) {
  allforgamma<-c()
  for(i in (1:1)){
    x<-c(rgamma(5400, shape=a/100, rate = 1))
    targetgam<-lmrgam(para = c(a/100, 1), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    mean1<-mean(x)
    tm1<-mean(x,1/9)
    wm1<-mean(winsor(x,fraction=1/9))
    sd1<-sd(x) 
    rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
    rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    #rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=FALSE,subsample=54000,sorted=TRUE)
    rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    targetm<-targetgam[1]
    targetvar<-(a/100)
    targettm<-((sqrt(a/100))^3)*2/sqrt(a/100)
    targetfm<-((sqrt(a/100))^4)*((6/(a/100))+3)
    targetl2<-targetgam[2]
    targetl3<-targetgam[3]*targetgam[2]
    targetl4<-targetgam[4]*targetgam[2]
    all<-c(
      meanBias=((mean1-targetm)/(sd1)),
      tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1[3]-targetm)/(sd1)),rmBias=((rqmmm1[2]-targetm)/(sd1)),qmBias=((rqmmm1[3]-targetm)/(sd1)),
      rvarBias=((rqscale1[4]-targetvar)/(rqscale1[6])),qvarBias=((rqscale1[5]-targetvar)/(rqscale1[6])),rl2Bias=((rqscale1[1]-targetl2)/(rqscale1[3])),ql2Bias=((rqscale1[2]-targetl2)/(rqscale1[3])),
      rtmBias=((rqtm1[4]-targettm)/(rqtm1[6])),qtmBias=((rqtm1[5]-targettm)/(rqtm1[6])),rl3Bias=((rqtm1[1]-targetl3)/(rqtm1[3])),ql3Bias=((rqtm1[2]-targetl3)/(rqtm1[3])),
      rfmBias=(rqfm1[4]-targetfm)/(rqfm1[6]),qfmBias=(rqfm1[5]-targetfm)/(rqfm1[6]),rl4Bias=((rqfm1[1]-targetl4)/(rqfm1[3])),ql4Bias=((rqfm1[2]-targetl4)/(rqfm1[3])))
    allforgamma<-rbind(allforgamma,all)
  }
  allforgamma[is.infinite(allforgamma)] <-NA
  listgamma<-rbind(listgamma,allforgamma)
}

write.csv(listgamma,paste("Congamma",batchnumber,".csv", sep = ","), row.names = TRUE)

#better for heavy-tailed distributions, because the kurtosis much larger than 4
listlnorm<-data.frame()
for (a in (1:300)) {
  allforlnorm<-c()
  for(i in (1:1)){
    x<-c(rlnorm(5400,meanlog=0,sdlog=a/100))
    targetlnorm<-lmrln3(para = c(0,0, a/100), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    mean1<-mean(x)
    tm1<-mean(x,1/9)
    wm1<-mean(winsor(x,fraction=1/9))
    sd1<-sd(x) 
    rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
    rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    #rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=FALSE,subsample=54000,sorted=TRUE)
    rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    targetm<-targetlnorm[1]
    targetvar<-(exp((a/100)^2)*(-1+exp((a/100)^2)))
    targettm<-sqrt(exp((a/100)^2)-1)*((2+exp((a/100)^2)))*((sqrt(exp((a/100)^2)*(-1+exp((a/100)^2))))^3)
    targetfm<-((-3+exp(4*((a/100)^2))+2*exp(3*((a/100)^2))+3*exp(2*((a/100)^2))))*((sqrt(exp((a/100)^2)*(-1+exp((a/100)^2))))^4)
    targetl2<-targetlnorm[2]
    targetl3<-targetlnorm[3]*targetlnorm[2]
    targetl4<-targetlnorm[4]*targetlnorm[2]
    all<-c(
      meanBias=((mean1-targetm)/(sd1)),
      tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1[3]-targetm)/(sd1)),rmBias=((rqmmm1[2]-targetm)/(sd1)),qmBias=((rqmmm1[3]-targetm)/(sd1)),
      rvarBias=((rqscale1[4]-targetvar)/(rqscale1[6])),qvarBias=((rqscale1[5]-targetvar)/(rqscale1[6])),rl2Bias=((rqscale1[1]-targetl2)/(rqscale1[3])),ql2Bias=((rqscale1[2]-targetl2)/(rqscale1[3])),
      rtmBias=((rqtm1[4]-targettm)/(rqtm1[6])),qtmBias=((rqtm1[5]-targettm)/(rqtm1[6])),rl3Bias=((rqtm1[1]-targetl3)/(rqtm1[3])),ql3Bias=((rqtm1[2]-targetl3)/(rqtm1[3])),
      rfmBias=(rqfm1[4]-targetfm)/(rqfm1[6]),qfmBias=(rqfm1[5]-targetfm)/(rqfm1[6]),rl4Bias=((rqfm1[1]-targetl4)/(rqfm1[3])),ql4Bias=((rqfm1[2]-targetl4)/(rqfm1[3])))
    allforlnorm<-rbind(allforlnorm,all)
  }
  allforlnorm[is.infinite(allforlnorm)] <-NA
  listlnorm<-rbind(listlnorm,allforlnorm)
}

write.csv(listlnorm,paste("Conlnorm",batchnumber,".csv", sep = ","), row.names = TRUE)

listpareto<-data.frame()
for (a in (1:300)) {
  allforpareto<-c()
  for(i in (1:1)){
    library(VGAM)
    x<-c(VGAM::rpareto(5400, scale  = 1, shape=2+a/100))
    targetlpareto<-lmrgpa(para = c(1,1/(2+a/100),- 1/(2+a/100)), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    mean1<-mean(x)
    tm1<-mean(x,1/9)
    wm1<-mean(winsor(x,fraction=1/9))
    sd1<-sd(x) 
    rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
    rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    #rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=FALSE,subsample=54000,sorted=TRUE)
    rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    targetm<-targetlpareto[1]
    targetvar<-(((2+a/100))*(1)/((-2+(2+a/100))*((-1+(2+a/100))^2)))
    targettm<-((((2+a/100)+1)*(2)*(sqrt(a/100)))/((-3+(2+a/100))*(((2+a/100))^(1/2))))*(((sqrt(((2+a/100))*(1)/((-2+(2+a/100))*((-1+(2+a/100))^2))))^3))
    targetfm<-(3+(6*((2+a/100)^3+(2+a/100)^2-6*(2+a/100)-2)/(((2+a/100))*((-3+(2+a/100)))*((-4+(2+a/100))))))*((sqrt(((2+a/100))*(1)/((-2+(2+a/100))*((-1+(2+a/100))^2))))^4)
    targetl2<-targetlpareto[2]
    targetl3<-targetlpareto[3]*targetlpareto[2]
    targetl4<-targetlpareto[4]*targetlpareto[2]
    all<-c(
      meanBias=((mean1-targetm)/(sd1)),
      tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1[3]-targetm)/(sd1)),rmBias=((rqmmm1[2]-targetm)/(sd1)),qmBias=((rqmmm1[3]-targetm)/(sd1)),
      rvarBias=((rqscale1[4]-targetvar)/(rqscale1[6])),qvarBias=((rqscale1[5]-targetvar)/(rqscale1[6])),rl2Bias=((rqscale1[1]-targetl2)/(rqscale1[3])),ql2Bias=((rqscale1[2]-targetl2)/(rqscale1[3])),
      rtmBias=((rqtm1[4]-targettm)/(rqtm1[6])),qtmBias=((rqtm1[5]-targettm)/(rqtm1[6])),rl3Bias=((rqtm1[1]-targetl3)/(rqtm1[3])),ql3Bias=((rqtm1[2]-targetl3)/(rqtm1[3])),
      rfmBias=(rqfm1[4]-targetfm)/(rqfm1[6]),qfmBias=(rqfm1[5]-targetfm)/(rqfm1[6]),rl4Bias=((rqfm1[1]-targetl4)/(rqfm1[3])),ql4Bias=((rqfm1[2]-targetl4)/(rqfm1[3])))
    allforpareto<-rbind(allforpareto,all)
  }
  allforpareto[is.infinite(allforpareto)] <-NA
  listpareto<-rbind(listpareto,allforpareto)
}

write.csv(listpareto,paste("Conpareto",batchnumber,".csv", sep = ","), row.names = TRUE)
