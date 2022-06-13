

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

finddtm<-function (x,interval=9,fast=TRUE,batch=1000,sorted=TRUE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  subtract<-t(combn(sortedx, 3))
  getm<-function(vector){ 
    ((1/6)*(2*vector[1]-vector[2]-vector[3])*(-1*vector[1]+2*vector[2]-vector[3])*(-vector[1]-vector[2]+2*vector[3]))
  }
  dp2m<-apply(subtract,MARGIN=1,FUN=getm)
  getlm<-function(vector){ 
    (vector[3]-2*vector[2]+vector[1])
  }
  dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
  dps<-sort(dp2lm,decreasing = FALSE,method ="radix")
  dp2s<-sort(dp2m,decreasing = FALSE,method ="radix")
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
    dqtm<-log(((quatileexpectdp2s-mx1dp2s)/(abs(1-mx1dp2s)*((mx1dp2s-mx2dp2s)*2))),base=((abs(mx1dp2s-mx2dp2s)*2)))
  }else{
    dqtm<-log(((quatileexpectdp2s-mx1dp2s)/(abs(0-mx1dp2s)*((mx1dp2s-mx2dp2s)*2))),base=((abs(mx1dp2s-mx2dp2s)*2)))
  }
  drtm<-(expectdp2s-etm1dp2s[2])/(etm1dp2s[2]-etm1dp2s[3])
  if (mx1dps>0.5){
    dql3<-log(((quatileexpectdps-mx1dps)/(abs(1-mx1dps)*((mx1dps-mx2dps)*2))),base=((abs(mx1dps-mx2dps)*2)))
  }else{
    dql3<-log(((quatileexpectdps-mx1dps)/(abs(0-mx1dps)*((mx1dps-mx2dps)*2))),base=((abs(mx1dps-mx2dps)*2)))
  }
  drl3<-(expectdps-etm1dps[2])/(etm1dps[2]-etm1dps[3])
  all<-c(dqtm,drtm,dql3,drl3)
  return(all)
}
finddfm<-function (x,interval=9,fast=TRUE,batch=1000,sorted=TRUE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  subtract<-t(combn(sortedx, 4))
  getm<-function(vector){ 
    resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                  4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                  vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                  6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                  vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                               6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
    (resd)
  }
  dp2m<-apply(subtract,MARGIN=1,FUN=getm)
  getlm<-function(vector){ 
    (vector[4]-3*vector[3]+3*vector[2]-vector[1])
  }
  dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
  dps<-sort(dp2lm,decreasing = FALSE,method ="radix")
  dp2s<-sort(dp2m,decreasing = FALSE,method ="radix")
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
  
  mx1dp2s
  if (mx1dp2s>0.5){
    dqfm<-log(((quatileexpectdp2s-mx1dp2s)/(abs(1-mx1dp2s)*((mx1dp2s-mx2dp2s)*2))),base=((abs(mx1dp2s-mx2dp2s)*2)))
  }else{
    dqfm<-log(((quatileexpectdp2s-mx1dp2s)/(abs(0-mx1dp2s)*((mx1dp2s-mx2dp2s)*2))),base=((abs(mx1dp2s-mx2dp2s)*2)))
  }
  drfm<-(expectdp2s-etm1dp2s[2])/(etm1dp2s[2]-etm1dp2s[3])
  if (mx1dps>0.5){
    dql4<-log(((quatileexpectdps-mx1dps)/(abs(1-mx1dps)*((mx1dps-mx2dps)*2))),base=((abs(mx1dps-mx2dps)*2)))
  }else{
    dql4<-log(((quatileexpectdps-mx1dps)/(abs(0-mx1dps)*((mx1dps-mx2dps)*2))),base=((abs(mx1dps-mx2dps)*2)))
  }
  drl4<-(expectdps-etm1dps[2])/(etm1dps[2]-etm1dps[3])
  all<-c(dqfm,drfm,dql4,drl4)
  return(all)
}
dataslicing<-function (x,slicesize=90,FUN=finddtm){
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

#notice that quantile l3 moments for exponential distribution is not valid since the true mean is within etm and median

simulatedbatch<-c()
for(i in (1:10)){
  x<-c(rexp(5400,1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  dtm1<-dataslicing(x,slicesize=216,FUN=finddtm)
  dfm1<-dataslicing(x,slicesize=90,FUN=finddfm)
  all<-c(dqtm=dtm1[1],drtm=dtm1[2],dql3=dtm1[3],drl3=dtm1[4],dqfm=dfm1[1],drfm=dfm1[2],dql4=dfm1[3],drl4=dfm1[4])
  simulatedbatch<-rbind(simulatedbatch,all)
}
simulatedbatch[is.infinite(simulatedbatch)] <-NA

write.csv(simulatedbatch,paste("simulatedd",batchnumber,".csv", sep = ","), row.names = TRUE)

