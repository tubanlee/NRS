
greatest_common_divisor<- function(a, b) {
  if (b == 0) a else Recall(b, a %% b)
}
least_common_multiple<-function(a,b){
  g<-greatest_common_divisor(a, b)
  return(a/g * b)
}
data_augmentation<-function (x,targetsize){
  lengthx<-length(x)
  meanx<-mean(x)
  disn<-least_common_multiple(lengthx,targetsize)
  orderedx<-sort(x,decreasing = FALSE,method ="radix")
  disori<-rep(orderedx, each = disn/lengthx)
  group<-rep(1:targetsize, each =disn/targetsize)
  data_augmentationresult<-sapply(split(disori, group), mean)
  return(data_augmentationresult)
}

#equinveral trimmed mean and complementary trimmed mean
etm<-function (x,interval=9,fast=TRUE,batch="auto"){
  lengthx<-length(x)
  if (batch=="auto" ){
    batch<-ceiling(500000/lengthx)
  }
  Ksamples<-lengthx/interval
  IntKsamples<-ceiling(Ksamples)
  target1<-IntKsamples*interval
  if (Ksamples%%1!=0 ){
    if (fast==TRUE & lengthx<10000){
      x_ordered<-data_augmentation(x,target1)
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
      return((c((mean(batch1)),mean(batch2),mean(batch3))))
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
    return(Groupmean)}
}
rlaplace<-function (n,location,scale) {
  sample1<-runif(n)
  sample1<-location - sign(sample1 - 0.5) * scale * (log(2) + ifelse(sample1 < 0.5, log(sample1), log1p(-sample1)))
  sample1[scale <= 0] <- NaN
  sample1
}
rRayleigh<-function (n, scale) {
  sample1 <- scale * sqrt(-2 * log(runif(n)))
  sample1[scale <= 0] <- NaN
  sample1
}
rpareto<-function (n, scale, shape) {
  sample1 <- scale*(runif(n))^(-1/shape)
  sample1[scale <= 0] <- NaN
  sample1[shape <= 0] <- NaN
  sample1
}
mmme<-function(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  etm1<-etm(sortedx,interval=interval,fast=fast,batch=batch)
  mx1<-(min(which(sortedx>(etm1[2])))-1)/length(x)
  mx2<-1/2
  if (mx1>0.5){
    quatiletarget<-abs(1-mx1)*((mx1-mx2)*2)*(((abs(mx1-mx2)*2))^dqm)+mx1
  }else{
    quatiletarget<-abs(0-mx1)*((mx1-mx2)*2)*(((abs(mx1-mx2)*2))^dqm)+mx1
  }
  upper1<-(1-1/interval)
  lower1<-1/interval
  if (!is.na(quatiletarget) & quatiletarget>(upper1)){
    print(paste("Warning: the percentile exceeds ",as.character(upper1*interval),"/",as.character(interval),", the robustness shrinks."))
  }else if(!is.na(quatiletarget) & quatiletarget<(lower1)){
    print(paste("Warning: the percentile exceeds ",as.character(lower1*interval),"/",as.character(interval),", the robustness shrinks."))
  }
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  qm1<-quantile(sortedx,quatiletarget)
  rm1<--drm*etm1[3]+etm1[2]+drm*etm1[2]
  names(rm1)<-NULL
  listd<-c(mean(sortedx),etm1[2],rm1,qm1,tm1,wm1)
  return(listd)
}

rqmean<-function (x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = FALSE,alpha = 0.05,nboot = 1000){
  if(cise){
    return (mmmeci(x, interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm,alpha = alpha,nboot = nboot))
  } else {return (mmme(x,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm))}
}

mmmeci<-function(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,alpha=0.05,nboot=1000){
  data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  pairs1 <- data.frame()
  pairs2 <- data.frame()
  pairs3 <- data.frame()
  pairs4 <- data.frame()
  pairs5 <- data.frame()
  pairs6 <- data.frame()
  for (i in 1:nboot) {
    robustlocation1<-mmme(data[i,],interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
    pairs1 <- rbind(pairs1,robustlocation1[1])
    pairs2 <- rbind(pairs2,robustlocation1[2])
    pairs3 <- rbind(pairs3,robustlocation1[3])
    pairs4 <- rbind(pairs4,robustlocation1[4])
    pairs5 <- rbind(pairs3,robustlocation1[5])
    pairs6 <- rbind(pairs4,robustlocation1[6])
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
  estimate=mmme(x,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  estimate=c(mean=estimate[1],etm=estimate[2],rm=estimate[3],qm=estimate[4],tm1=estimate[5],wm1=estimate[6])
  result <- list(cimean=c(bootlist1[low],bootlist1[up]), cietm=c(bootlist2[low],bootlist2[up]),cirm=c(bootlist3[low],bootlist3[up]),
                 ciqm=c(bootlist4[low],bootlist4[up]),citm=c(bootlist5[low],bootlist5[up]),
                 ciwm=c(bootlist6[low],bootlist6[up]),semean=sd(bootlist1),seetm=sd(bootlist2),serm=sd(bootlist3),seqm=sd(bootlist4),setm=sd(bootlist5),sewm=sd(bootlist6),estimate=estimate)
  return(result)
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
for(i in (1:100)){
  x<-c(rlaplace(5400, location = 0, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-rqmean(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = TRUE,alpha = 0.05,nboot = 1000)
  targetm<-0
  all<-c(cimean=c(rqmmm1$cimean),cietm=c(rqmmm1$cietm),cirm=c(rqmmm1$cirm),ciqm=c(rqmmm1$ciqm),citm=c(rqmmm1$citm),ciwm=c(rqmmm1$ciwm),semean=rqmmm1$semean,seetm=rqmmm1$seetm,serm=rqmmm1$serm,seqm=rqmmm1$seqm,setm=rqmmm1$setm,sewm=rqmmm1$sewm,mean=rqmmm1$estimate[1],
         etm=rqmmm1$estimate[2],rm=rqmmm1$estimate[3],qm=rqmmm1$estimate[4],tm=rqmmm1$estimate[5],wm=rqmmm1$estimate[6],sd=sd1,meanBias=((rqmmm1$estimate[1]-targetm)/(sd1)),
    tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1$estimate[2]-targetm)/(sd1)),rmBias=((rqmmm1$estimate[3]-targetm)/(sd1)),qmBias=((rqmmm1$estimate[4]-targetm)/(sd1))
    )
  allforlaplace<-rbind(allforlaplace,all)
}
allforlaplace[is.infinite(allforlaplace)] <-NA
colMeans(allforlaplace)
write.csv(allforlaplace,paste("Conlaplace,",batchnumber,".csv", sep = ","), row.names = TRUE)


#the biases of robust/quantile fourth moments for normal are huge, because the kurtosis is 3, too differ from that of exponential distribution (9)

#but robust 4th moment still remain high consistent, because the distributions of U-statistics of L-moments are more symmetric.

allfornorm<-c()
for(i in (1:100)){
  x<-c(rnorm(5400))
  x<-sort(x,decreasing = FALSE,method ="radix")
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-rqmean(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = TRUE,alpha = 0.05,nboot = 1000)
  targetm<-0
  all<-c(cimean=c(rqmmm1$cimean),cietm=c(rqmmm1$cietm),cirm=c(rqmmm1$cirm),ciqm=c(rqmmm1$ciqm),citm=c(rqmmm1$citm),ciwm=c(rqmmm1$ciwm),semean=rqmmm1$semean,seetm=rqmmm1$seetm,serm=rqmmm1$serm,seqm=rqmmm1$seqm,setm=rqmmm1$setm,sewm=rqmmm1$sewm,mean=rqmmm1$estimate[1],
         etm=rqmmm1$estimate[2],rm=rqmmm1$estimate[3],qm=rqmmm1$estimate[4],tm=rqmmm1$estimate[5],wm=rqmmm1$estimate[6],sd=sd1,meanBias=((rqmmm1$estimate[1]-targetm)/(sd1)),
         tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1$estimate[2]-targetm)/(sd1)),rmBias=((rqmmm1$estimate[3]-targetm)/(sd1)),qmBias=((rqmmm1$estimate[4]-targetm)/(sd1))
  )
  allfornorm<-rbind(allfornorm,all)
}
allfornorm[is.infinite(allfornorm)] <-NA
colMeans(allfornorm)
write.csv(allfornorm,paste("Connorm,",batchnumber,".csv", sep = ","), row.names = TRUE)

#the biases of robust/quantile fourth moments for logis are not very large ~0.2 and ~0.06, because the kurtosis is 4.2, better than normal

allforlogis<-c()
for(i in (1:100)){
  x<-c(rlogis(5400, location = 0, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-rqmean(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = TRUE,alpha = 0.05,nboot = 1000)
  targetm<-0
  all<-c(cimean=c(rqmmm1$cimean),cietm=c(rqmmm1$cietm),cirm=c(rqmmm1$cirm),ciqm=c(rqmmm1$ciqm),citm=c(rqmmm1$citm),ciwm=c(rqmmm1$ciwm),semean=rqmmm1$semean,seetm=rqmmm1$seetm,serm=rqmmm1$serm,seqm=rqmmm1$seqm,setm=rqmmm1$setm,sewm=rqmmm1$sewm,mean=rqmmm1$estimate[1],
         etm=rqmmm1$estimate[2],rm=rqmmm1$estimate[3],qm=rqmmm1$estimate[4],tm=rqmmm1$estimate[5],wm=rqmmm1$estimate[6],sd=sd1,meanBias=((rqmmm1$estimate[1]-targetm)/(sd1)),
         tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1$estimate[2]-targetm)/(sd1)),rmBias=((rqmmm1$estimate[3]-targetm)/(sd1)),qmBias=((rqmmm1$estimate[4]-targetm)/(sd1))
  )
  allforlogis<-rbind(allforlogis,all)
}
allforlogis[is.infinite(allforlogis)] <-NA
colMeans(allforlogis)
write.csv(allforlogis,paste("Conlogis,",batchnumber,".csv", sep = ","), row.names = TRUE)

#the biases of robust/quantile fourth moments for Rayleigh are very large ~0.3 and ~0.11, because the kurtosis is 3.245, very close to normal

allforRayleigh<-c()
for(i in (1:100)){
  x<-c(rrayleigh(5400, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-rqmean(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = TRUE,alpha = 0.05,nboot = 1000)
  targetm<-sqrt(pi/2)
  all<-c(cimean=c(rqmmm1$cimean),cietm=c(rqmmm1$cietm),cirm=c(rqmmm1$cirm),ciqm=c(rqmmm1$ciqm),citm=c(rqmmm1$citm),ciwm=c(rqmmm1$ciwm),semean=rqmmm1$semean,seetm=rqmmm1$seetm,serm=rqmmm1$serm,seqm=rqmmm1$seqm,setm=rqmmm1$setm,sewm=rqmmm1$sewm,mean=rqmmm1$estimate[1],
         etm=rqmmm1$estimate[2],rm=rqmmm1$estimate[3],qm=rqmmm1$estimate[4],tm=rqmmm1$estimate[5],wm=rqmmm1$estimate[6],sd=sd1,meanBias=((rqmmm1$estimate[1]-targetm)/(sd1)),
         tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1$estimate[2]-targetm)/(sd1)),rmBias=((rqmmm1$estimate[3]-targetm)/(sd1)),qmBias=((rqmmm1$estimate[4]-targetm)/(sd1))
  )
  allforRayleigh<-rbind(allforRayleigh,all)
}
allforRayleigh[is.infinite(allforRayleigh)] <-NA
colMeans(allforRayleigh)
write.csv(allforRayleigh,paste("ConRayleigh,",batchnumber,".csv", sep = ","), row.names = TRUE)


allforexp<-c()
for(i in (1:100)){
  x<-c(rexp(5400,1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  tm1<-mean(x,1/9)
  wm1<-mean(winsor(x,fraction=1/9))
  sd1<-sd(x) 
  rqmmm1<-rqmean(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = TRUE,alpha = 0.05,nboot = 1000)
  targetm<-1
  all<-c(cimean=c(rqmmm1$cimean),cietm=c(rqmmm1$cietm),cirm=c(rqmmm1$cirm),ciqm=c(rqmmm1$ciqm),citm=c(rqmmm1$citm),ciwm=c(rqmmm1$ciwm),semean=rqmmm1$semean,seetm=rqmmm1$seetm,serm=rqmmm1$serm,seqm=rqmmm1$seqm,setm=rqmmm1$setm,sewm=rqmmm1$sewm,mean=rqmmm1$estimate[1],
         etm=rqmmm1$estimate[2],rm=rqmmm1$estimate[3],qm=rqmmm1$estimate[4],tm=rqmmm1$estimate[5],wm=rqmmm1$estimate[6],sd=sd1,meanBias=((rqmmm1$estimate[1]-targetm)/(sd1)),
         tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1$estimate[2]-targetm)/(sd1)),rmBias=((rqmmm1$estimate[3]-targetm)/(sd1)),qmBias=((rqmmm1$estimate[4]-targetm)/(sd1))
  )
  allforexp<-rbind(allforexp,all)
}
allforexp[is.infinite(allforexp)] <-NA
colMeans(allforexp)
write.csv(allforexp,paste("Conexp,",batchnumber,".csv", sep = ","), row.names = TRUE)


#same issue, when kurtosis lower than 4, the performances of robust/quantile fourth moments are very poor.


#the 0.01-0.1 shape parameter range were removed, as the variance too high, often introducing bugs.
library(lmom)
listWeibull<-data.frame()
for (a in (10:300)) {
  allforweibull<-c()
  for(i in (1:30)){
    x<-c(rweibull(5400, shape=a/100, scale = 1))
    targetwei<-lmrwei(para = c(0, 1, a/100), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    tm1<-mean(x,1/9)
    wm1<-mean(winsor(x,fraction=1/9))
    sd1<-sd(x) 
    rqmmm1<-rqmean(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = TRUE,alpha = 0.05,nboot = 1000)
    targetm<-gamma(1+1/(a/100))
    all<-c(cimean=c(rqmmm1$cimean),cietm=c(rqmmm1$cietm),cirm=c(rqmmm1$cirm),ciqm=c(rqmmm1$ciqm),citm=c(rqmmm1$citm),ciwm=c(rqmmm1$ciwm),semean=rqmmm1$semean,seetm=rqmmm1$seetm,serm=rqmmm1$serm,seqm=rqmmm1$seqm,setm=rqmmm1$setm,sewm=rqmmm1$sewm,mean=rqmmm1$estimate[1],
           etm=rqmmm1$estimate[2],rm=rqmmm1$estimate[3],qm=rqmmm1$estimate[4],tm=rqmmm1$estimate[5],wm=rqmmm1$estimate[6],sd=sd1,meanBias=((rqmmm1$estimate[1]-targetm)/(sd1)),
           tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1$estimate[2]-targetm)/(sd1)),rmBias=((rqmmm1$estimate[3]-targetm)/(sd1)),qmBias=((rqmmm1$estimate[4]-targetm)/(sd1))
    )
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
  for(i in (1:30)){
    x<-c(rgamma(5400, shape=a/100, rate = 1))
    targetgam<-lmrgam(para = c(a/100, 1), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    tm1<-mean(x,1/9)
    wm1<-mean(winsor(x,fraction=1/9))
    sd1<-sd(x) 
    rqmmm1<-rqmean(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = TRUE,alpha = 0.05,nboot = 1000)
    targetm<-targetgam[1]
    all<-c(cimean=c(rqmmm1$cimean),cietm=c(rqmmm1$cietm),cirm=c(rqmmm1$cirm),ciqm=c(rqmmm1$ciqm),citm=c(rqmmm1$citm),ciwm=c(rqmmm1$ciwm),semean=rqmmm1$semean,seetm=rqmmm1$seetm,serm=rqmmm1$serm,seqm=rqmmm1$seqm,setm=rqmmm1$setm,sewm=rqmmm1$sewm,mean=rqmmm1$estimate[1],
           etm=rqmmm1$estimate[2],rm=rqmmm1$estimate[3],qm=rqmmm1$estimate[4],tm=rqmmm1$estimate[5],wm=rqmmm1$estimate[6],sd=sd1,meanBias=((rqmmm1$estimate[1]-targetm)/(sd1)),
           tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1$estimate[2]-targetm)/(sd1)),rmBias=((rqmmm1$estimate[3]-targetm)/(sd1)),qmBias=((rqmmm1$estimate[4]-targetm)/(sd1))
    )
    
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
  for(i in (1:30)){
    x<-c(rlnorm(5400,meanlog=0,sdlog=a/100))
    targetlnorm<-lmrln3(para = c(0,0, a/100), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    tm1<-mean(x,1/9)
    wm1<-mean(winsor(x,fraction=1/9))
    sd1<-sd(x) 
    rqmmm1<-rqmean(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = TRUE,alpha = 0.05,nboot = 1000)
    targetm<-targetlnorm[1]
    all<-c(cimean=c(rqmmm1$cimean),cietm=c(rqmmm1$cietm),cirm=c(rqmmm1$cirm),ciqm=c(rqmmm1$ciqm),citm=c(rqmmm1$citm),ciwm=c(rqmmm1$ciwm),semean=rqmmm1$semean,seetm=rqmmm1$seetm,serm=rqmmm1$serm,seqm=rqmmm1$seqm,setm=rqmmm1$setm,sewm=rqmmm1$sewm,mean=rqmmm1$estimate[1],
           etm=rqmmm1$estimate[2],rm=rqmmm1$estimate[3],qm=rqmmm1$estimate[4],tm=rqmmm1$estimate[5],wm=rqmmm1$estimate[6],sd=sd1,meanBias=((rqmmm1$estimate[1]-targetm)/(sd1)),
           tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1$estimate[2]-targetm)/(sd1)),rmBias=((rqmmm1$estimate[3]-targetm)/(sd1)),qmBias=((rqmmm1$estimate[4]-targetm)/(sd1))
    )
    allforlnorm<-rbind(allforlnorm,all)
  }
  allforlnorm[is.infinite(allforlnorm)] <-NA
  listlnorm<-rbind(listlnorm,allforlnorm)
}

write.csv(listlnorm,paste("Conlnorm",batchnumber,".csv", sep = ","), row.names = TRUE)

listpareto<-data.frame()
for (a in (1:300)) {
  allforpareto<-c()
  for(i in (1:30)){
    x<-c(rpareto(5400, scale  = 1, shape=2+a/100))
    targetlpareto<-lmrgpa(para = c(1,1/(2+a/100),- 1/(2+a/100)), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    tm1<-mean(x,1/9)
    wm1<-mean(winsor(x,fraction=1/9))
    sd1<-sd(x) 
    rqmmm1<-rqmean(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = TRUE,alpha = 0.05,nboot = 1000)
    targetm<-targetlpareto[1]
    all<-c(cimean=c(rqmmm1$cimean),cietm=c(rqmmm1$cietm),cirm=c(rqmmm1$cirm),ciqm=c(rqmmm1$ciqm),citm=c(rqmmm1$citm),ciwm=c(rqmmm1$ciwm),semean=rqmmm1$semean,seetm=rqmmm1$seetm,serm=rqmmm1$serm,seqm=rqmmm1$seqm,setm=rqmmm1$setm,sewm=rqmmm1$sewm,mean=rqmmm1$estimate[1],
           etm=rqmmm1$estimate[2],rm=rqmmm1$estimate[3],qm=rqmmm1$estimate[4],tm=rqmmm1$estimate[5],wm=rqmmm1$estimate[6],sd=sd1,meanBias=((rqmmm1$estimate[1]-targetm)/(sd1)),
           tmBias=((tm1-targetm)/(sd1)),wmBias=((wm1-targetm)/(sd1)), etmBias=((rqmmm1$estimate[2]-targetm)/(sd1)),rmBias=((rqmmm1$estimate[3]-targetm)/(sd1)),qmBias=((rqmmm1$estimate[4]-targetm)/(sd1))
    )
    allforpareto<-rbind(allforpareto,all)
  }
  allforpareto[is.infinite(allforpareto)] <-NA
  listpareto<-rbind(listpareto,allforpareto)
}

write.csv(listpareto,paste("Conpareto",batchnumber,".csv", sep = ","), row.names = TRUE)