




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
etm<-function (x,interval=9,fast=TRUE,batch="auto"){
  lengthx<-length(x)
  if (batch=="auto" ){
    batch<-ceiling(500000/lengthx)+1
  }
  Ksamples<-lengthx/interval
  IntKsamples<-ceiling(Ksamples)
  target1<-IntKsamples*interval
  etmass<-function(x1,IntKsamples,target1,sorted=FALSE){
    if(sorted){
      x_ordered<-x1
    }else{
      x_ordered<-sort(x1,decreasing = FALSE,method ="radix")
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
    return(Groupmean)
  }
  if (Ksamples%%1!=0 ){
    if (fast==TRUE & lengthx<10000){
      x_ordered<-data_augmentation(x,target1)
      Groupmean<-etmass(x_ordered,IntKsamples,target1,sorted=TRUE)
      return(Groupmean)
    }
    else{
      addedx<-matrix(sample(x,size=(target1-lengthx)*batch,replace=TRUE),nrow=batch)
      xt<-t(as.data.frame(x))
      xmatrix<-as.data.frame(lapply(xt, rep, batch))
      allmatrix<-cbind(addedx,xmatrix)
      batchresults<-apply(allmatrix,1,etmass,IntKsamples=IntKsamples,target1=target1,sorted=FALSE)
      return(colMeans(t(batchresults)))}
  }
  else{
    Groupmean<-etmass(x,IntKsamples,target1,sorted=FALSE)
    return(Groupmean)}
}

#load finite sample bias corrected d
fd <- read.csv(("resultsfd_Rayleigh.csv"))

finited<-function(n,fd,type){
  indextypelist<-c("rm","qm","rl2","ql2","rvar","qvar","rl3","ql3","rtm","qtm","rl4","ql4","rfm","qfm")
  indextype<-which(indextypelist==(type))+1
  if (type=="qm"){
    if(n%in% fd[,1]){
      return(fd[which(fd[,1] == n),indextype])
    }
    else if(n>5400){
      return(0.821497)
    }
    else{
      maxn<-max(which(fd[,1] < n))
      minn<-min(which(fd[,1] > n))
      size1<-fd[maxn,1]
      size2<-fd[minn,1]
      d1<-fd[maxn,indextype]
      d2<-fd[minn,indextype]
      return(((d2-d1)*((n-size1)/(size2-size1)))+d1)}}
  if (type=="rm"){
    if(n%in% fd[,1]){
      return(fd[which(fd[,1] == n),indextype])
    }
    else if(n>5400){
      return(0.366919)
    }
    else{
      maxn<-max(which(fd[,1] < n))
      minn<-min(which(fd[,1] > n))
      size1<-fd[maxn,1]
      size2<-fd[minn,1]
      d1<-fd[maxn,indextype]
      d2<-fd[minn,indextype]
      return(((d2-d1)*((n-size1)/(size2-size1)))+d1)}}
  else{
    if(n%in% fd[,1]){
      return(fd[which(fd[,1] == n),indextype])
    }
    else if(n>5400){
      return(fd[117,indextype])
    }
    else{
      maxn<-max(which(fd[,1] < n))
      minn<-min(which(fd[,1] > n))
      size1<-fd[maxn,1]
      size2<-fd[minn,1]
      d1<-fd[maxn,indextype]
      d2<-fd[minn,indextype]
      return(((d2-d1)*((n-size1)/(size2-size1)))+d1)}}
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
  drm<-finited(n=length(x),fd=fd,type="rm")
  dqm<-finited(n=length(x),fd=fd,type="qm")
  
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
  drm<-finited(n=length(x),fd=fd,type="rl2")
  dqm<-finited(n=length(x),fd=fd,type="ql2")
  dlmo<-mmm(expectboot=mean(dp2lm),expecttrue=expectdps,x=dp2lm,interval=9,fast=TRUE,batch=10000,sorted=FALSE,drm=drm,dqm=dqm)
  drm<-finited(n=length(x),fd=fd,type="rvar")
  dqm<-finited(n=length(x),fd=fd,type="qvar")
  dmo<-mmm(expectboot=mean(dp2m),expecttrue=expectdp2s,x=dp2m,interval=9,fast=TRUE,batch=10000,sorted=FALSE,drm=drm,dqm=dqm)
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
  drm<-finited(n=length(x),fd=fd,type="rl3")
  dqm<-finited(n=length(x),fd=fd,type="ql3")
  dlmo<-mmm(expectboot=mean(dp2lm),expecttrue=expectdps,x=dp2lm,interval=9,fast=fast,batch=batch,sorted=FALSE,drm=drm,dqm=dqm)
  drm<-finited(n=length(x),fd=fd,type="rtm")
  dqm<-finited(n=length(x),fd=fd,type="qtm")
  dmo<-mmm(expectboot=mean(dp2m),expecttrue=expectdp2s,x=dp2m,interval=9,fast=fast,batch=batch,sorted=FALSE,drm=drm,dqm=dqm)
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
  drm<-finited(n=lengthn,fd=fd,type="rl4")
  dqm<-finited(n=lengthn,fd=fd,type="ql4")
  dlmo<-mmm(expectboot=mean(dp2lm),expecttrue=expectdps,x=dp2lm,interval=9,fast=fast,batch=batch,sorted=FALSE,drm=drm,dqm=dqm)
  drm<-finited(n=lengthn,fd=fd,type="rfm")
  dqm<-finited(n=lengthn,fd=fd,type="qfm")
  dmo<-mmm(expectboot=mean(dp2m),expecttrue=expectdp2s,x=dp2m,interval=9,fast=fast,batch=batch,sorted=FALSE,drm=drm,dqm=dqm)
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
eLaplace<-function (n,location,scale) {
  sample1<-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))
  sample1<-location - sign(sample1 - 0.5) * scale * (log(2) + ifelse(sample1 < 0.5, log(sample1), log1p(-sample1)))
  sample1
}

library(Lmoments)
library(foreach)
library(doParallel)
numCores <- detectCores()
registerDoParallel(numCores) 
#parallel


#this finite sample consistency is just a qualitative evaluation, since equinterval simulation is biased.

allforlaplace<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
  library(VGAM)
  library(Lmoments)
  x<-c(eLaplace(n=i, location = 0, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  mean1<-mean(x)
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
  rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
  rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
  rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
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
}
allforlaplace[is.infinite(allforlaplace)] <-NA

write.csv(allforlaplace,paste("Conlaplace,",batchnumber,".csv", sep = ","), row.names = TRUE)


enorm<-function (n, location,scale) {
  library(pracma)
  sample1 <- location+(scale)*sqrt(2)*erfinv(2*(1-seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))-1)
  sample1
}
allfornorm<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
  library(Lmoments)
  x<-c(enorm(n=i,location-0,scale=1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  mean1<-mean(x)
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-mmme(x,interval=9,fast=TRUE,batch=1000,sorted=FALSE,drm=0.366439,dqm=0.8241268)
  rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
  rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
  rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
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
}
allfornorm[is.infinite(allfornorm)] <-NA

write.csv(allfornorm,paste("Connorm,",batchnumber,".csv", sep = ","), row.names = TRUE)

elogis<-function (n,location,scale) {
  sample1<-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))
  sample1<-location + scale * log((1-sample1)/sample1)
  sample1
}

allforlogis<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
  library(Lmoments)
  x<-c(elogis(n=i, location = 0, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  mean1<-mean(x)
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
  rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
  rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
  rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
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
}
allforlogis[is.infinite(allforlogis)] <-NA

write.csv(allforlogis,paste("Conlogis,",batchnumber,".csv", sep = ","), row.names = TRUE)

eRayleigh<-function (n, scale) {
  sample1 <- scale * sqrt(-2 * log((seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))))
  sample1[scale <= 0] <- NaN
  sample1
}
allforRayleigh<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
  library(VGAM)
  library(Lmoments)
  x<-c(eRayleigh(n=i, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  mean1<-mean(x)
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
  rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
  rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
  rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
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
}
allforRayleigh[is.infinite(allforRayleigh)] <-NA

write.csv(allforRayleigh,paste("ConRayleigh,",batchnumber,".csv", sep = ","), row.names = TRUE)

eexp<-function (n, scale) {
  sample1 <- (-1/scale)*(log(1-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))))
  sample1[scale <= 0] <- NaN
  sample1
}

allforexp<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
  library(Lmoments)
  x<-c(eexp(n=i,scale=1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  mean1<-mean(x)
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
  rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
  rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
  rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
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
}
allforexp[is.infinite(allforexp)] <-NA

write.csv(allforexp,paste("Conexp,",batchnumber,".csv", sep = ","), row.names = TRUE)


eWeibull<-function (n,shape, scale = 1){
  sample1<-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))
  sample1<-qweibull(sample1,shape=shape, scale = scale)
  sample1
}
library(lmom)
listWeibull<-data.frame()
for (a in (10:300)) {
  allforweibull<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
    library(lmom)
    library(Lmoments)
    x<-c(eWeibull(n=i, shape=a/100, scale = 1))
    targetwei<-lmrwei(para = c(0, 1, a/100), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    mean1<-mean(x)
    tm1<-mean(x,1/9)
    wm1<-mean(winsor(x,fraction=1/9))
    sd1<-sd(x) 
    rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
    rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
    rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
    rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
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
  }
  allforweibull[is.infinite(allforweibull)] <-NA
  listWeibull<-rbind(listWeibull,allforweibull)
}

write.csv(listWeibull,paste("ConWeibull",batchnumber,".csv", sep = ","), row.names = TRUE)
egamma<-function (n,shape,scale = 1) {
  sample1<-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))
  sample1<-qgamma(sample1,shape=shape,scale=scale)
  sample1
}

library(lmom)
listgamma<-data.frame()
for (a in (1:300)) {
  allforgamma<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
    library(lmom)
    library(Lmoments)
    x<-c(egamma(n=i, shape=a/100, rate = 1))
    targetgam<-lmrgam(para = c(a/100, 1), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    mean1<-mean(x)
    tm1<-mean(x,1/9)
    wm1<-mean(winsor(x,fraction=1/9))
    sd1<-sd(x) 
    rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
    rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
    rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
    rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
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
  }
  allforgamma[is.infinite(allforgamma)] <-NA
  listgamma<-rbind(listgamma,allforgamma)
}

write.csv(listgamma,paste("Congamma",batchnumber,".csv", sep = ","), row.names = TRUE)
elnorm<-function (n,location,scale) {
  library(pracma)
  sample1 <- exp(location+sqrt(scale^2)*sqrt(2)*erfinv(2*(1-seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))-1))
  sample1
}
listlnorm<-data.frame()
for (a in (1:300)) {
  allforlnorm<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
    library(lmom)
    library(Lmoments)
    x<-c(elnorm(n=i,location=0,scale=a/100))
    targetlnorm<-lmrln3(para = c(0,0, a/100), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    mean1<-mean(x)
    tm1<-mean(x,1/9)
    wm1<-mean(winsor(x,fraction=1/9))
    sd1<-sd(x) 
    rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
    rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
    rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
    rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
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
  }
  allforlnorm[is.infinite(allforlnorm)] <-NA
  listlnorm<-rbind(listlnorm,allforlnorm)
}

write.csv(listlnorm,paste("Conlnorm",batchnumber,".csv", sep = ","), row.names = TRUE)
ePareto<-function (n, scale, shape) {
  sample1 <- scale*((seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1))))^(-1/shape)
  sample1[scale <= 0] <- NaN
  sample1[shape <= 0] <- NaN
  sample1
}
listpareto<-data.frame()
for (a in (1:300)) {
  allforpareto<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
    library(lmom)
    library(VGAM)
    library(Lmoments)
    x<-c(ePareto(n=i, scale  = 1, shape=2+a/100))
    targetlpareto<-lmrgpa(para = c(1,1/(2+a/100),- 1/(2+a/100)), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    mean1<-mean(x)
    tm1<-mean(x,1/9)
    wm1<-mean(winsor(x,fraction=1/9))
    sd1<-sd(x) 
    rqmmm1<-mmm(expectboot=mean1,expecttrue=mean1,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
    rqscale1<-rqscale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
    rqtm1<-rqtm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
    rqfm1<-rqfm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=5400,sorted=TRUE)
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
  }
  allforpareto[is.infinite(allforpareto)] <-NA
  listpareto<-rbind(listpareto,allforpareto)
}

write.csv(listpareto,paste("Conpareto",batchnumber,".csv", sep = ","), row.names = TRUE)