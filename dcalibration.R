




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
finddmmm<-function(x,interval=9,fast=TRUE,batch=10000,sorted=TRUE){
  expect<-mean(x)
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  quatileexpect<-(min(which(sortedx>(expect)))-1)/length(x)
  etm1<-etm(x,interval=interval,fast=fast,batch=batch)
  mx1<-(min(which(sortedx>(etm1[2])))-1)/length(x)
  mx2<-1/2
  if (mx1>0.5){
    qm1<-log(((quatileexpect-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
  }else{
    qm1<-log(((quatileexpect-mx1)/(abs(0-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
  }
  rm1<-(expect-etm1[2])/(etm1[2]-etm1[3])
  listd<-c(qm1,rm1)
  return(listd)
}
finddscale<-function (x,interval=9,fast=TRUE,batch=1000,sorted=TRUE){
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
  dp2s<-sort(dp2,decreasing = FALSE,method ="radix")
  etm1dps<-etm(dps,interval=interval,fast=fast,batch=batch)
  etm1dp2s<-etm(dp2s,interval=interval,fast=fast,batch=batch)
  expectdps<-mean(dps)
  expectdp2s<-mean(dp2s)
  quatileexpectdps<-(min(which(dps>(expectdps)))-1)/length(dps)
  mx1dps<-(min(which(dps>(etm1dps[2])))-1)/length(dps)
  mx2dps<-1/2
  quatileexpectdp2s<-(min(which(dp2s>(expectdp2s)))-1)/length(dp2s)
  mx1dp2s<-(min(which(dp2s>(etm1dp2s[2])))-1)/length(dp2s)
  mx2dp2s<-1/2
  if (mx1dp2s>0.5){
    dqsd<-log(((quatileexpectdp2s-mx1dp2s)/(abs(1-mx1dp2s)*((mx1dp2s-mx2dp2s)*2))),base=((abs(mx1dp2s-mx2dp2s)*2)))
  }else{
    dqsd<-log(((quatileexpectdp2s-mx1dp2s)/(abs(0-mx1dp2s)*((mx1dp2s-mx2dp2s)*2))),base=((abs(mx1dp2s-mx2dp2s)*2)))
  }
  drsd<-(expectdp2s-etm1dp2s[2])/(etm1dp2s[2]-etm1dp2s[3])
  if (mx1dps>0.5){
    dql2<-log(((quatileexpectdps-mx1dps)/(abs(1-mx1dps)*((mx1dps-mx2dps)*2))),base=((abs(mx1dps-mx2dps)*2)))
  }else{
    dql2<-log(((quatileexpectdps-mx1dps)/(abs(0-mx1dps)*((mx1dps-mx2dps)*2))),base=((abs(mx1dps-mx2dps)*2)))
  }
  drl2<-(expectdps-etm1dps[2])/(etm1dps[2]-etm1dps[3])
  all<-c(dqsd,drsd,dql2,drl2)
  return(all)
}

simulatedbatch<-c()
for(i in (1:10)){
  x<-c(rexp(5400,1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  dmmm1<-finddmmm(x,interval=9,fast=TRUE,batch=10000,sorted=TRUE)
  dscale1<-finddscale(x,interval=9,fast=TRUE,batch=1000,sorted=TRUE)
  all<-c(dqm=dmmm1[1],drm=dmmm1[2],dqsd=dscale1[1],drsd=dscale1[2],dql2=dscale1[3],drl2=dscale1[4])
  simulatedbatch<-rbind(simulatedbatch,all)
}
simulatedbatch[is.infinite(simulatedbatch)] <-NA

write.csv(simulatedbatch,paste("simulatedd",batchnumber,".csv", sep = ","), row.names = TRUE)
