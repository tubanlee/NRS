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
#load finite sample bias corrected d
fd <- read.csv(("resultsfd.csv"))
finited<-function(n,fd,type){
  if (type=="qm"){
    if(n%in% fd$size){
      return(fd$dqm[which(fd$size == n)])
    }
    else if(n>5400){
      return(0.821497)
    }
    else{
      maxn<-max(which(fd$size < n))
      minn<-min(which(fd$size > n))
      size1<-fd$size[maxn]
      size2<-fd$size[minn]
      d1<-fd$dqm[maxn]
      d2<-fd$dqm[minn]
      return(((d2-d1)*((n-size1)/(size2-size1)))+d1)}}
  if (type=="rm"){
    if(n%in% fd$size){
      return(fd$drm[which(fd$size == n)])
    }else if(n>5400){
      return(0.366919)
    }
    else{
      maxn<-max(which(fd$size < n))
      minn<-min(which(fd$size > n))
      size1<-fd$size[maxn]
      size2<-fd$size[minn]
      d1<-fd$drm[maxn]
      d2<-fd$drm[minn]
      return(((d2-d1)*((n-size1)/(size2-size1)))+d1)}}
  if (type=="qsd"){
    if(n%in% fd$size){
      return(fd$dqsd[which(fd$size == n)])
    }else if(n>5400){
      return(fd$dqsd[which(fd$size == 5400)])
    }
    else{
      maxn<-max(which(fd$size < n))
      minn<-min(which(fd$size > n))
      size1<-fd$size[maxn]
      size2<-fd$size[minn]
      d1<-fd$dqsd[maxn]
      d2<-fd$dqsd[minn]
      return(((d2-d1)*((n-size1)/(size2-size1)))+d1)}}
  if (type=="rsd"){
    if(n%in% fd$size){
      return(fd$drsd[which(fd$size == n)])
    }else if(n>5400){
      return(fd$drsd[which(fd$size == 5400)])
    }
    
    else{
      maxn<-max(which(fd$size < n))
      minn<-min(which(fd$size > n))
      size1<-fd$size[maxn]
      size2<-fd$size[minn]
      d1<-fd$drsd[maxn]
      d2<-fd$drsd[minn]
      return(((d2-d1)*((n-size1)/(size2-size1)))+d1)}}
  if (type=="ql2"){
    if(n%in% fd$size){
      return(fd$dql2[which(fd$size == n)])
    }else if(n>5400){
      return(fd$dql2[which(fd$size == 5400)])
    }
    else{
      maxn<-max(which(fd$size < n))
      minn<-min(which(fd$size > n))
      size1<-fd$size[maxn]
      size2<-fd$size[minn]
      d1<-fd$dql2[maxn]
      d2<-fd$dql2[minn]
      return(((d2-d1)*((n-size1)/(size2-size1)))+d1)}}
  if (type=="rl2"){
    if(n%in% fd$size){
      return(fd$drl2[which(fd$size == n)])
    }else if(n>5400){
      return(fd$drl2[which(fd$size == 5400)])
    }else{
      maxn<-max(which(fd$size < n))
      minn<-min(which(fd$size > n))
      size1<-fd$size[maxn]
      size2<-fd$size[minn]
      d1<-fd$drl2[maxn]
      d2<-fd$drl2[minn]
      return(((d2-d1)*((n-size1)/(size2-size1)))+d1)}}
}
#quantile mean
qm<-function(x,interval=9,fast=TRUE,batch=1000,d=0.821497,sorted=FALSE,fsbc=FALSE){
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  samplesize<-length(x)
  etmx<-c(etm(sortedx,interval=interval,fast=fast,batch=batch))
  mx1<-(min(which(sortedx>(etmx[2])))-1)/length(x)
  mx2<-1/2
  if (fsbc){
    d<-finited(n=samplesize,fd=fd,type="qm")
  }
  quatiletarget<-abs(1-mx1)*((mx1-mx2)*2)*(((abs(mx1-mx2)*2))^d)+mx1
  if (quatiletarget>8/9){
    print("Warning: the percentile exceeds 8/9, the robustness shrinks")
  }else if(quatiletarget<1/9){
    print("Warning: the percentile exceeds 1/9, the robustness shrinks")
  }
  resultdx<-quantile(sortedx,quatiletarget)
  (resultdx)
}
#robust mean
rm<-function(x,interval=9,fast=TRUE,batch=10000,d=0.366919,fsbc=FALSE){
  samplesize<-length(x)
  if (fsbc){
    d<-finited(n=samplesize,fd=fd,type="rm")
  }
  etm1<-etm(x,interval=interval,fast=fast,batch=batch)
  -d*etm1[3]+etm1[2]+d*etm1[2] 
}
#quantile standard deviation
qsd<-function (x,interval=9,fast=FALSE,batch=540000,d=0.78609163,sorted=FALSE,fsbc=FALSE){
  if (length(x)>5000& fast==FALSE){
    print("Warning: The computational complexity is (e*n/2)^2, data slicing or bootstrap is recommended")
  }
  samplesize<-length(x)
  if (fsbc){
    d<-finited(n=samplesize,fd=fd,type="qsd")
  }
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  get<-function(vector){ 
    ((vector[1]-vector[2])^2)/2
  }
  if (fast){
    subtract<-t(replicate(batch, sample(sortedx, size = 2)))
    dp2<-apply(subtract,MARGIN=1,FUN=get)
  }else{
    subtract<-sapply(sortedx, "-", sortedx)
    subtract[lower.tri(subtract)] <- NA
    diag(subtract)=NA
    subtract<-na.omit(as.vector(subtract))
    dp<-subtract[subtract>0]
    dp2<-(dp^2)/2
  }
  resultdp<-(qm(dp2,interval=interval,fast=TRUE,batch=batch,d=d))
  (sqrt(resultdp))
}
#robust standard deviation
rsd<-function (x,interval=9,fast=FALSE,batch=540000,d=0.794006982,sorted=FALSE,fsbc=FALSE){
  if (length(x)>5000& fast==FALSE){
    print("Warning: The computational complexity is (e*n/2)^2, data slicing or bootstrap is recommended")
  }
  samplesize<-length(x)
  if (fsbc){
    d<-finited(n=samplesize,fd=fd,type="rsd")
  }
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  get<-function(vector){ 
    ((vector[1]-vector[2])^2)/2
  }
  if (fast){
    subtract<-t(replicate(batch, sample(sortedx, size = 2)))
    dp2<-apply(subtract,MARGIN=1,FUN=get)
  }else{
    subtract<-sapply(sortedx, "-", sortedx)
    subtract[lower.tri(subtract)] <- NA
    diag(subtract)=NA
    subtract<-na.omit(as.vector(subtract))
    dp<-subtract[subtract>0]
    dp2<-(dp^2)/2
  }
  resultdp<-(rm(dp2,interval=interval,fast=TRUE,batch=batch,d=d))
  sqrt(resultdp)
}
#quantile L2-moment
ql2<-function (x,interval=9,fast=FALSE,batch=540000,d=0.824464716,sorted=FALSE,fsbc=FALSE){
  if (length(x)>5000& fast==FALSE){
    print("Warning: The computational complexity is (e*n/2)^2, data slicing or bootstrap is recommended")
  }
  samplesize<-length(x)
  if (fsbc){
    d<-finited(n=samplesize,fd=fd,type="ql2")
  }
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  get<-function(vector){ 
    abs(vector[2]-vector[1])/2
  }
  if (fast){
    subtract<-t(replicate(batch, sample(sortedx, size = 2)))
    dp<-apply(subtract,MARGIN=1,FUN=get)
  }else{
    subtract<-sapply(sortedx, "-", sortedx)
    subtract[lower.tri(subtract)] <- NA
    diag(subtract)=NA
    subtract<-na.omit(as.vector(subtract))
    dp<-(subtract[subtract>0])/2
  }
  resultdp<-(qm(dp,interval=interval,fast=TRUE,batch=batch,d=d))
  (resultdp)
}
#robust L2-moment
rl2<-function (x,interval=9,fast=FALSE,batch=540000,d=0.366885946,sorted=FALSE,fsbc=FALSE){
  if (length(x)>5000& fast==FALSE){
    print("Warning: The computational complexity is n^2, data slicing or bootstrap is recommended")
  }
  samplesize<-length(x)
  if (fsbc){
    d<-finited(n=samplesize,fd=fd,type="rl2")
  }
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  get<-function(vector){ 
    abs(vector[1]-vector[2])/2
  }
  if (fast){
    subtract<-t(replicate(batch, sample(sortedx, size = 2)))
    dp<-apply(subtract,MARGIN=1,FUN=get)
  }else{
    subtract<-sapply(sortedx, "-", sortedx)
    subtract[lower.tri(subtract)] <- NA
    diag(subtract)=NA
    subtract<-na.omit(as.vector(subtract))
    dp<-(subtract[subtract>0])/2
  }
  resultdp<-(rm(dp,interval=interval,fast=TRUE,batch=batch,d=d))
  (resultdp)
}

#quantile skewness
qskew<-function (x,interval=9,fast=FALSE,batch=540000,d=0.58778,sorted=FALSE){
  if (length(x)>300& fast==FALSE){
    print("Warning: The computational complexity is n^3, data slicing or bootstrap is recommended")
  }
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  get<-function(vector){ 
    ((1/6)*(2*vector[1]-vector[2]-vector[3])*(-1*vector[1]+2*vector[2]-vector[3])*(-vector[1]-vector[2]+2*vector[3]))
  }
  if (fast){
    subtract<-t(replicate(batch, sort(sample(sortedx, size = 3))))
  }else{
    subtract<-t(combn(sortedx, 3))
  }
  dp2<-apply(subtract,MARGIN=1,FUN=get)
  resultdp<-(qm(dp2,interval=interval,fast=TRUE,batch=batch,d=d))
  sd1<-qsd(x,interval=9,fast=TRUE,fsbc=TRUE,batch=540000)
  (resultdp)/((sd1)^3)
}
#robust skewness
rskew<-function (x,interval=9,fast=TRUE,batch=540000,d=1.71172,sorted=FALSE){
  if (length(x)>300& fast==FALSE){
    print("Warning: The computational complexity is n^3, data slicing or bootstrap is recommended")
  }
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  get<-function(vector){ 
    ((1/6)*(2*vector[1]-vector[2]-vector[3])*(-1*vector[1]+2*vector[2]-vector[3])*(-vector[1]-vector[2]+2*vector[3]))
  }
  if (fast){
    subtract<-t(replicate(batch, sort(sample(sortedx, size = 3))))
  }else{
    subtract<-t(combn(sortedx, 3))
  }
  dp2<-apply(subtract,MARGIN=1,FUN=get)
  resultdp<-(rm(dp2,interval=interval,fast=TRUE,batch=batch,d=d))
  sd1<-rsd(x,interval=9,fast=TRUE,fsbc=TRUE,batch=540000)
  (resultdp)/((sd1[1])^3)
}
#quantile scaled L3-moment
ql3<-function (x,interval=9,fast=TRUE,batch=540000,d=1.18983,sorted=FALSE){
  if (length(x)>300& fast==FALSE){
    print("Warning: The computational complexity is n^3, data slicing or bootstrap is recommended")
  }
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  get<-function(vector){ 
    (vector[3]-2*vector[2]+vector[1])
  }
  if (fast){
    subtract<-t(replicate(batch, sort(sample(sortedx, size = 3))))
  }else{
    subtract<-t(combn(sortedx, 3))
  }
  dp2<-apply(subtract,MARGIN=1,FUN=get)
  resultdp<-(qm(dp2,interval=interval,fast=TRUE,batch=batch,d=d))
  l21<-ql2(x,interval=9,fast=TRUE,fsbc=TRUE,batch=540000)
  (((1/3)*resultdp)/((l21[1])))
}
#robust scaled L3-moment
rl3<-function (x,interval=9,fast=TRUE,batch=540000,d=0.176434,sorted=FALSE){
  if (length(x)>300& fast==FALSE){
    print("Warning: The computational complexity is n^3, data slicing or bootstrap is recommended")
  }
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  get<-function(vector){ 
    (vector[3]-2*vector[2]+vector[1])
  }
  if (fast){
    subtract<-t(replicate(batch, sort(sample(sortedx, size = 3))))
  }else{
    subtract<-t(combn(sortedx, 3))
  }
  dp2<-apply(subtract,MARGIN=1,FUN=get)
  resultdp<-(rm(dp2,interval=interval,fast=TRUE,batch=batch,d=d))
  l21<-rl2(x,interval=9,fast=TRUE,fsbc=TRUE,batch=540000)
  (((1/3)*resultdp)/((l21[1])))
}
#quantile scaled L4-moment
ql4<-function (x,interval=9,fast=TRUE,batch=540000,d=1.511272749 ,sorted=FALSE){
  if (length(x)>108& fast==FALSE){
    print("Warning: The computational complexity is n^4, data slicing or bootstrap is recommended")
  }
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  get<-function(vector){ 
    (vector[4]-3*vector[3]+3*vector[2]-vector[1])
  }
  if (fast){
    subtract<-t(replicate(batch, sort(sample(sortedx, size = 4))))
  }else{
    subtract<-t(combn(sortedx, 4))
  }
  dp2<-apply(subtract,MARGIN=1,FUN=get)
  resultdp<-(qm(dp2,interval=interval,fast=TRUE,batch=batch,d=d))
  l21<-ql2(x,interval=9,fast=TRUE,fsbc=TRUE,batch=540000)
  (((1/4)*resultdp)/((l21[1])))
}
#robust scaled L4-moment
rl4<-function (x,interval=9,fast=TRUE,batch=540000,d=0.009055351,sorted=FALSE){
  if (length(x)>108& fast==FALSE){
    print("Warning: The computational complexity is n^4, data slicing or bootstrap is recommended")
  }
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  get<-function(vector){ 
    (vector[4]-3*vector[3]+3*vector[2]-vector[1])
  }
  if (fast){
    subtract<-t(replicate(batch, sort(sample(sortedx, size = 4))))
  }else{
    subtract<-t(combn(sortedx, 4))
  }
  dp2<-apply(subtract,MARGIN=1,FUN=get)
  resultdp<-(rm(dp2,interval=interval,fast=TRUE,batch=batch,d=d))
  l21<-rl2(x,interval=9,fast=TRUE,fsbc=TRUE,batch=540000)
  (((1/4)*resultdp)/((l21[1])))
}
#quantile kurtosis
qkurt<-function (x,interval=9,fast=TRUE,batch=540000,d=0.1364722 ,sorted=FALSE){
  if (length(x)>108& fast==FALSE){
    print("Warning: The computational complexity is n^4, data slicing or bootstrap is recommended")
  }
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  get<-function(vector){ 
    resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                  4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                  vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                  6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                  vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                               6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
    (resd)
  }
  if (fast){
    subtract<-t(replicate(batch, sort(sample(sortedx, size = 4))))
  }else{
    subtract<-t(combn(sortedx, 4))
  }
  dp2<-apply(subtract,MARGIN=1,FUN=get)
  resultdp<-(qm(dp2,interval=interval,fast=TRUE,batch=batch,d=d))
  sd1<-qsd(x,interval=9,fast=TRUE,fsbc=TRUE,batch=540000)
  (qkurt=(resultdp)/((sd1[1])^4))
}
#robust kurtosis
rkurt<-function (x,interval=9,fast=TRUE,batch=540000,d=3.378723,sorted=FALSE){
  if (length(x)>108& fast==FALSE){
    print("Warning: The computational complexity is n^4, data slicing or bootstrap is recommended")
  }
  if (sorted){
    sortedx<-x
  }else {
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  get<-function(vector){ 
    resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                  4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                  vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                  6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                  vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                               6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
    (resd)
  }
  if (fast){
    subtract<-t(replicate(batch, sort(sample(sortedx, size = 4))))
  }else{
    subtract<-t(combn(sortedx, 4))
  }
  dp2<-apply(subtract,MARGIN=1,FUN=get)
  resultdp<-(rm(dp2,interval=interval,fast=TRUE,batch=batch,d=d))
  sd1<-rsd(x,interval=9,fast=TRUE,fsbc=TRUE,batch=540000)
  (resultdp)/((sd1[1])^4)
}
#reducing time complexity to nlogn with dataslicing
dataslicing<-function (x,slicesize=108,FUN=rskew){
  lengthx<-length(x)
  slices<-lengthx/slicesize
  IntKslices<-floor(slices)
  if (IntKslices==1){
    return("sample size too small, no need for data slicing")
  }
  x_ordered<-sort(x,decreasing = FALSE,method ="radix")
  group<-rep(rep(c(1:IntKslices), each=1), times=slicesize)
  suppressWarnings(Group1<-split(x_ordered, group))
  Groupall<-((sapply(Group1,FUN)))
  if (is.null(dim(Groupall))){
    mean(Groupall)
  }else{
    apply(Groupall,MARGIN=1,FUN=mean)}
}
#reducing time complexity to nlogn with dataslicing
dataslicing2<-function (x,slicesize=108,FUN=rskew,interval=9,fast=TRUE,batch=1000,d=NULL,sorted=FALSE){
  if (is.null(d)){
    return("please provide d")
  }
  lengthx<-length(x)
  slices<-lengthx/slicesize
  IntKslices<-floor(slices)
  x_ordered<-sort(x,decreasing = FALSE,method ="radix")
  group<-rep(rep(c(1:IntKslices), each=1), times=slicesize)
  suppressWarnings(Group1<-split(x_ordered, group))
  Groupall0<-c()
  for (i in (1:IntKslices)){
    Groupall<-FUN(unlist(Group1[i]),interval=interval,fast=fast,batch=batch,d=d,sorted=sorted)
    Groupall0<-rbind(Groupall0,Groupall)
  }
  if (is.null(dim(Groupall0))){
    all1<-mean(Groupall0)
  }else{
    all1<-apply(Groupall0,MARGIN=2,FUN=mean)}
  (all1)
}

#test
x<-rexp(5400,1)
#the population mean is 1
mean(x)
rm(x,interval=9,fast=TRUE,fsbc=FALSE)
#finite sample bias correction
rm(x,interval=9,fast=TRUE,fsbc=TRUE)

qm(x,interval=9,fast=TRUE,fsbc=FALSE)
qm(x,interval=9,fast=TRUE,fsbc=TRUE)
#the population standard deviation is 1
rsd(x,interval=9,fast=FALSE,fsbc=TRUE)
rsd(x,interval=9,fast=TRUE,fsbc=TRUE,batch=540000)
qsd(x,interval=9,fast=FALSE,fsbc=TRUE)
#bootstrap 540000 for around three decimal accuracy, the batch should be the factor of interval
qsd(x,interval=9,fast=TRUE,fsbc=TRUE,batch=540000)

#the population L2-moment is 1/2
rl2(x,interval=9,fast=FALSE,fsbc=TRUE)
rl2(x,interval=9,fast=TRUE,fsbc=TRUE,batch=540000)
ql2(x,interval=9,fast=FALSE,fsbc=TRUE)
ql2(x,interval=9,fast=TRUE,fsbc=TRUE,batch=540000)

#the population skewness is 2
rskew(x,interval=9,fast=TRUE,batch=540000)
qskew(x,interval=9,fast=TRUE,batch=540000)
#the population L-skewness is 1/3
rl3(x,interval=9,fast=TRUE,batch=540000)
ql3(x,interval=9,fast=TRUE,batch=540000)

#the population kurtosis is 9
rkurt(x,interval=9,fast=TRUE,batch=540000)
qkurt(x,interval=9,fast=TRUE,batch=540000)
#the population L-kurtosis is 1/6
rl4(x,interval=9,fast=TRUE,batch=540000)
ql4(x,interval=9,fast=TRUE,batch=540000)

#test paramater defined dataslicing
dataslicing2(x,slicesize=90,FUN=rskew,interval=9,fast=TRUE,batch=1000,d=1.71172,sorted=FALSE)

x<-rgamma(5400,4,1)
#the population mean is 4
mean(x)
rm(x)
qm(x)
#the population standard deviation is 2
rsd(x)
qsd(x)
#the population L2-moment is 1.09375
rl2(x)
ql2(x)
#the population skewness is 1
rskew(x,interval=9,fast=TRUE,batch=540000)
qskew(x,interval=9,fast=TRUE,batch=540000)
#the population L-skewness is 0.1646599 
rl3(x,interval=9,fast=TRUE,batch=540000)
ql3(x,interval=9,fast=TRUE,batch=540000)

#the population kurtosis is 4.5
rkurt(x,interval=9,fast=TRUE,batch=540000)
qkurt(x,interval=9,fast=TRUE,batch=540000)
#the population L-kurtosis is 0.1312521 
rl4(x,interval=9,fast=TRUE,batch=540000)
ql4(x,interval=9,fast=TRUE,batch=540000)

x<-rpois(5400,8)
#the population mean is 8
mean(x)
rm(x)
#quantile mean returns 7 or 9
qm(x)
#the population standard deviation is sqrt(8)
rsd(x)
qsd(x)
#the population L2-moment is 1.583
rl2(x)
#quantile l2 returns 2
ql2(x)
#the population skewness is 0.3535534
rskew(x,interval=9,fast=TRUE,batch=540000)
qskew(x,interval=9,fast=TRUE,batch=540000)
#the population L-skewness is 0.0592
rl3(x,interval=9,fast=TRUE,batch=540000)
ql3(x,interval=9,fast=TRUE,batch=540000)

#the population kurtosis is 1/8+3
rkurt(x,interval=9,fast=TRUE,batch=540000)
qkurt(x,interval=9,fast=TRUE,batch=540000)
#the population L-kurtosis is 0.1204 
rl4(x,interval=9,fast=TRUE,batch=540000)
ql4(x,interval=9,fast=TRUE,batch=540000)

#quantile L-moments is not suitable for discrete unimodal distributions

#for more tests, use the codes in consistency.R
