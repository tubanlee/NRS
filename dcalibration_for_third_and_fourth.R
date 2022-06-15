



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
finddmmm<-function(expect,x,interval=9,fast=TRUE,batch=10000,sorted=FALSE){
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

finddtm<-function (x,interval=9,fast=TRUE,dataslicing=FALSE,batch=540000,sorted=TRUE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  if (fast){
    subtract<-t(replicate(batch, sort(sample(sortedx, size = 3))))
  }else if(dataslicing){
    subtract<-sortedx
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
  
  if (dataslicing){
    dp2lm<-dataslicingm(x=subtract,slicesize=18,FUN=alllm)
  }else{
    dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
  }
  
  getm<-function(vector){ 
    ((1/6)*(2*vector[1]-vector[2]-vector[3])*(-1*vector[1]+2*vector[2]-vector[3])*(-vector[1]-vector[2]+2*vector[3]))
  }
  allm<-function(sortedx){ 
    subtract<-t(combn(sortedx, 3))
    apply(subtract,MARGIN=1,FUN=getm)
  }
  if (dataslicing){
    dp2m<-dataslicingm(x=subtract,slicesize=18,FUN=allm)
  }else{
    dp2m<-apply(subtract,MARGIN=1,FUN=getm)
  }

  lengthn<-length(sortedx)
  lm1<-Lmoments(sortedx)
  expectdps<-lm1[3]
  expectdp2s<-((sum((sortedx - mean(sortedx))^3)/lengthn)*(lengthn^2/((lengthn-1)*(lengthn-2))))
  
  dlmo<-finddmmm(expect=expectdps,x=dp2lm,interval=9,fast=TRUE,batch=10000,sorted=FALSE)
  dmo<-finddmmm(expect=expectdp2s,x=dp2m,interval=9,fast=TRUE,batch=10000,sorted=FALSE)
  dqtm<-dmo[1]
  drtm<-dmo[2]
  dql3<-dlmo[1]
  drl3<-dlmo[2]
  all<-c(dqtm,drtm,dql3,drl3)
  return(all)
}
finddfm<-function (x,interval=9,fast=TRUE,dataslicing=FALSE,batch=1000,sorted=TRUE){
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
  if (fast){
    subtract<-t(replicate(batch, sort(sample(sortedx, size = 4))))
  }else if(dataslicing){
    subtract<-sortedx
  }else{
    subtract<-t(combn(sortedx, 4))
  }
  if (dataslicing){
    dp2m<-dataslicingm(x=subtract,slicesize=18,FUN=allm)
  }else{
    dp2m<-apply(subtract,MARGIN=1,FUN=getm)
  }
  getlm<-function(vector){ 
    (1/4)*(vector[4]-3*vector[3]+3*vector[2]-vector[1])
  }
  alllm<-function(sortedx){ 
    subtract<-t(combn(sortedx, 4))
    apply(subtract,MARGIN=1,FUN=getlm)
  }
  if (dataslicing){
    dp2lm<-dataslicingm(x=subtract,slicesize=18,FUN=alllm)
  }else{
    dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
  }

  lengthn<-length(sortedx)
  lm1<-Lmoments(sortedx)
  expectdps<-lm1[4]
  expectdp2s<-(sum((sortedx - mean(sortedx))^4)/lengthn)
  
  dlmo<-finddmmm(expect=expectdps,x=dp2lm,interval=9,fast=TRUE,batch=10000,sorted=FALSE)
  dmo<-finddmmm(expect=expectdp2s,x=dp2m,interval=9,fast=TRUE,batch=10000,sorted=FALSE)
  dqfm<-dmo[1]
  drfm<-dmo[2]
  dql4<-dlmo[1]
  drl4<-dlmo[2]
  all<-c(dqfm,drfm,dql4,drl4)
  return(all)
}
dataslicingm<-function (x,slicesize=108,FUN=rskew){
  lengthx<-length(x)
  slices<-lengthx/slicesize
  IntKslices<-floor(slices)
  x_ordered<-sort(x,decreasing = FALSE,method ="radix")
  group<-rep(rep(c(1:IntKslices), each=1), times=slicesize)
  suppressWarnings(Group1<-split(x_ordered, group))
  Groupall<-((sapply(Group1,FUN)))
  Groupall
}

simulatedbatch<-c()
for(i in (1:10)){
  x<-c(rexp(5400,1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  dtm1<-finddtm(x,interval=9,fast=TRUE,batch=540000,sorted=TRUE)
  dfm1<-finddfm(x,interval=9,fast=TRUE,batch=540000,sorted=TRUE)
  all<-c(dqtm=dtm1[1],drtm=dtm1[2],dql3=dtm1[3],drl3=dtm1[4],dqfm=dfm1[1],drfm=dfm1[2],dql4=dfm1[3],drl4=dfm1[4])
  simulatedbatch<-rbind(simulatedbatch,all)
}

simulatedbatch[is.infinite(simulatedbatch)] <-NA

write.csv(simulatedbatch,paste("simulatedd",batchnumber,".csv", sep = ","), row.names = TRUE)


simulatedbatchdataslicing<-c()
for(i in (1:10)){
  x<-c(rexp(5400,1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  dtm1<-finddtm(x,interval=9,fast=FALSE,batch=540000,sorted=TRUE,dataslicing=TRUE)
  dfm1<-finddfm(x,interval=9,fast=FALSE,batch=540000,sorted=TRUE,dataslicing=TRUE)
  all<-c(dqtm=dtm1[1],drtm=dtm1[2],dql3=dtm1[3],drl3=dtm1[4],dqfm=dfm1[1],drfm=dfm1[2],dql4=dfm1[3],drl4=dfm1[4])
  simulatedbatchdataslicing<-rbind(simulatedbatchdataslicing,all)
}

simulatedbatchdataslicing[is.infinite(simulatedbatchdataslicing)] <-NA

write.csv(simulatedbatchdataslicing,paste("simulatedddataslicing",batchnumber,".csv", sep = ","), row.names = TRUE)

