





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
#quantile mean
qm<-function(x,interval=9,fast=TRUE,batch=1000,d=0.821497,sorted=FALSE){
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  etmx<-c(etm(sortedx,interval=interval,fast=fast,batch=batch))
  mx1<-(min(which(sortedx>(etmx[2])))-1)/length(x)
  mx2<-1/2
  if (mx1>0.5){
    quatiletarget<-abs(1-mx1)*((mx1-mx2)*2)*(((abs(mx1-mx2)*2))^d)+mx1
  }else{
    quatiletarget<-abs(0-mx1)*((mx1-mx2)*2)*(((abs(mx1-mx2)*2))^d)+mx1
  }
  resultdx<-quantile(sortedx,quatiletarget)
  (resultdx)
}
#robust mean
rm<-function(x,interval=9,fast=TRUE,batch=10000,d=0.366919){
  etm1<-etm(x,interval=interval,fast=fast,batch=batch)
  -d*etm1[3]+etm1[2]+d*etm1[2] 
}
#standardized bias of robust/quantile standard deviation and L2 moment, all combined
Bias_rqscale<-function (x,target1,target2,interval=9,fast=TRUE,batch=1000,sorted=TRUE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  subtract<-sapply(sortedx, "-", sortedx)
  subtract[lower.tri(subtract)] <- NA
  diag(subtract)=NA
  subtract<-na.omit(as.vector(subtract))
  dp<-subtract[subtract>0]
  dp2<-(dp^2)/2
  dps<-sort(dp,decreasing = FALSE,method ="radix")
  dps<-dps/2
  dp2s<-sort(dp2,decreasing = FALSE,method ="radix")
  dp2s<-(dp2s)
  etm1dps<-etm(dps,interval=interval,fast=fast,batch=batch)
  etm1dp2s<-etm(dp2s,interval=interval,fast=fast,batch=batch)
  mx1dps<-(min(which(dps>(etm1dps[2])))-1)/length(dps)
  mx2dps<-1/2
  mx1dp2s<-(min(which(dp2s>(etm1dp2s[2])))-1)/length(dp2s)
  mx2dp2s<-1/2
  dqsd<-0.78609163
  dql2<-0.824464716
  drsd<-0.794006982
  drl2<-0.366885946
  if (mx1dp2s>0.5){
    quatiletarget2<-abs(1-mx1dp2s)*((mx1dp2s-mx2dp2s)*2)*(((abs(mx1dp2s-mx2dp2s)*2))^dqsd)+mx1dp2s
  }else{
    quatiletarget2<-abs(0-mx1dp2s)*((mx1dp2s-mx2dp2s)*2)*(((abs(mx1dp2s-mx2dp2s)*2))^dqsd)+mx1dp2s
  }
  if (mx1dps>0.5){
    quatiletarget<-abs(1-mx1dps)*((mx1dps-mx2dps)*2)*(((abs(mx1dps-mx2dps)*2))^dql2)+mx1dps
  }else{
    quatiletarget<-abs(0-mx1dps)*((mx1dps-mx2dps)*2)*(((abs(mx1dps-mx2dps)*2))^dql2)+mx1dps
  }
  qsd<-quantile(dp2s,quatiletarget2)
  rsd<--drsd*etm1dp2s[3]+etm1dp2s[2]+drsd*etm1dp2s[2]
  ql2<-quantile(dps,quatiletarget)
  rl2<--drl2*etm1dps[3]+etm1dps[2]+drl2*etm1dps[2] 
  standardizedbias<-c((sqrt(qsd)-target1)/sd(dp2s),(sqrt(rsd)-target1)/sd(dp2s),(ql2-target2)/sd(dps),(rl2-target2)/sd(dps))
  standardizedbias
}
winsor<-function (x, fraction=.05)
{
  lim <- quantile(x, probs=c(fraction, 1-fraction))
  x[ x < lim[1] ] <- lim[1]
  x[ x > lim[2] ] <- lim[2]
  x
}

allfornorm<-c()
for(i in (1:5)){
  x<-c(rnorm(5400,0))
  x<-sort(x,decreasing = FALSE,method ="radix")
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  qm1<-qm(x,interval=9,fast=TRUE,batch=1000,d=0.821497,sorted=TRUE)
  rm1<-rm(x,interval=9,fast=TRUE,batch=10000,d=0.366919)
  tm1Bias<-((((tm1))-0))/sd(x)  
  wm1Bias<-((((wm1))-0))/sd(x)
  qm1Bias<-((((qm1))-0))/sd(x)  
  rm1Bias<-((((rm1))-0))/sd(x)
  scale1<-Bias_rqscale(x,target1=1,target2=1/sqrt(pi),interval=9,fast=TRUE,batch=1000,sorted=TRUE)
  all<-c(tmBias=tm1Bias,wmBias=wm1Bias,qmBias=qm1Bias,rmBias=rm1Bias,sqsd=scale1[1],srsd=scale1[2,2],sql2=scale1[2,3],srl2=scale1[2,4])
  allfornorm<-rbind(allfornorm,all)
}
allfornorm[is.infinite(allfornorm)] <-NA

write.csv(allfornorm,paste("Connorm,",batchnumber,".csv", sep = ","), row.names = TRUE)


allforlaplace<-c()
for(i in (1:5)){
  library(VGAM)
  x<-c(rlaplace(5400, location = 0, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  qm1<-qm(x,interval=9,fast=TRUE,batch=1000,d=0.821497,sorted=TRUE)
  rm1<-rm(x,interval=9,fast=TRUE,batch=10000,d=0.366919)
  tm1Bias<-((((tm1))-0))/sd(x)  
  wm1Bias<-((((wm1))-0))/sd(x)
  qm1Bias<-((((qm1))-0))/sd(x)  
  rm1Bias<-((((rm1))-0))/sd(x)
  scale1<-Bias_rqscale(x,target1=sqrt(2),target2=3/4,interval=9,fast=TRUE,batch=1000,sorted=TRUE)
  all<-c(tmBias=tm1Bias,wmBias=wm1Bias,qmBias=qm1Bias,rmBias=rm1Bias,sqsd=scale1[1],srsd=scale1[2,2],sql2=scale1[2,3],srl2=scale1[2,4])
  allforlaplace<-rbind(allforlaplace,all)
}
allforlaplace[is.infinite(allforlaplace)] <-NA
write.csv(allforlaplace,paste("Conlaplace,",batchnumber,".csv", sep = ","), row.names = TRUE)

allforlogis<-c()
for(i in (1:5)){
  x<-c(rlogis(5400, location = 0, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  qm1<-qm(x,interval=9,fast=TRUE,batch=1000,d=0.821497,sorted=TRUE)
  rm1<-rm(x,interval=9,fast=TRUE,batch=10000,d=0.366919)
  tm1Bias<-((((tm1))-0))/sd(x)  
  wm1Bias<-((((wm1))-0))/sd(x)
  qm1Bias<-((((qm1))-0))/sd(x)  
  rm1Bias<-((((rm1))-0))/sd(x)
  scale1<-Bias_rqscale(x,target1=sqrt((pi^2)/3),target2=1,interval=9,fast=TRUE,batch=1000,sorted=TRUE)
  all<-c(tmBias=tm1Bias,wmBias=wm1Bias,qmBias=qm1Bias,rmBias=rm1Bias,sqsd=scale1[1],srsd=scale1[2,2],sql2=scale1[2,3],srl2=scale1[2,4])
  allforlogis<-rbind(allforlogis,all)
}
allforlogis[is.infinite(allforlogis)] <-NA

write.csv(allforlogis,paste("Conlogis,",batchnumber,".csv", sep = ","), row.names = TRUE)