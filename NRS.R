
#NRS

#I combined all the estimators into one function (easy for reviewing). There might be errors, and these are not bugs, 
#but because R is prone to producing errors for such a large function. 
#Run one function each time. If there is an error, try to restart, and then it will be fixed. 
#It will be completely fixed in the future by rewriting the code in C++.

#require library "Lmoments" to compare the consistencies and standard errors of NRSs.
if (!require("Lmoments")) install.packages("Lmoments")
library(Lmoments)

#require foreach and doparallel for parallel processing of bootstrap
if (!require("foreach")) install.packages("foreach")
library(foreach)
if (!require("doParallel")) install.packages("doParallel")
library(doParallel)
#registering clusters, can set a smaller number using numCores-1 
numCores <- detectCores()
registerDoParallel(numCores) 

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
mmm<-function(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224){
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
  qm1<-quantile(sortedx,quatiletarget)
  rm1<--drm*etm1[3]+etm1[2]+drm*etm1[2]
  names(rm1)<-NULL
  listd<-c(rm1,qm1)
  return(listd)
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
  qm1<-quantile(sortedx,quatiletarget)
  rm1<--drm*etm1[3]+etm1[2]+drm*etm1[2]
  names(rm1)<-NULL
  listd<-c(mean(sortedx),etm1[2],rm1,qm1)
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
  for (i in 1:nboot) {
    robustlocation1<-mmme(data[i,],interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
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
  estimate=mmme(x,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  estimate=c(mean=estimate[1],etm=estimate[2],rm=estimate[3],qm=estimate[4])
  result <- list(cimean=c(bootlist1[low],bootlist1[up]), cietm=c(bootlist2[low],bootlist2[up]),cirm=c(bootlist3[low],bootlist3[up]),
                 ciqm=c(bootlist4[low],bootlist4[up]),semean=sd(bootlist1),seetm=sd(bootlist2),serm=sd(bootlist3),seqm=sd(bootlist4),estimate=estimate)
  return(result)
}

rqscale<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.3659,drm=0.7930,dlqm=0.8218,dqm=0.7825){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  
  lengthn<-length(sortedx)
  if (boot){
    subtract<-t(replicate(times , sort(sample(sortedx, size = 2))))
    getlm<-function(vector){ 
      (vector[2]-vector[1])/2
    }
    dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
    getm<-function(vector){ 
      ((vector[1]-vector[2])^2)/2
    }
    dp2m<-apply(subtract,MARGIN=1,FUN=getm)
  }else{
    if (lengthn>5000){
      print("Warning: The computational complexity is (e*n/2)^2, bootstrap is recommended")
    }
    subtract<-sapply(sortedx, "-", sortedx)
    subtract[lower.tri(subtract)] <- NA
    diag(subtract)=NA
    subtract<-na.omit(as.vector(subtract))
    dp<-subtract[subtract>0]
    dp2lm<-dp/2
    dp2m<-(dp^2)/2
  }
  lm1<-Lmoments(sortedx)
  expectdps<-lm1[2]
  expectdp2s<-(sd(sortedx))^2
  dlmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
  dmo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  all<-c(expectdps,dlmo[1],dlmo[2],sd(dp2lm),
         expectdp2s,dmo[1],dmo[2],sd(dp2m)
  )
  return(all)
}

rqscaleci<-function(x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.3659,drm=0.7930,dlqm=0.8218,dqm=0.7825,alpha=0.05,nboot=1000){
  data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  pairs1 <- data.frame()
  pairs2 <- data.frame()
  pairs3 <- data.frame()
  pairs4 <- data.frame()
  pairs5 <- data.frame()
  pairs6 <- data.frame()
  for (i in 1:nboot) {
    robustlocation1<-rqscale(data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrm,drm=drm,dlqm=dlqm,dqm=dqm)
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
  estimate=rqscale(x,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrm,drm=drm,dlqm=dlqm,dqm=dqm)
  result <- list(cil2=c(bootlist1[low],bootlist1[up]), cirl2=c(bootlist2[low],bootlist2[up]),ciql2=c(bootlist3[low],bootlist3[up]),
                 civar=c(bootlist4[low],bootlist4[up]),cirvar=c(bootlist5[low],bootlist5[up]),ciqvar=c(bootlist6[low],bootlist6[up]),
                 sel2=sd(bootlist1),serl2=sd(bootlist2),seql2=sd(bootlist3),sevar=sd(bootlist4),servar=sd(bootlist5),seqvar=sd(bootlist6),
                 estimate=estimate)
  return(result)
}

rqtm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.1808,drm=1.7492,dlqm=1.1753,dqm=0.5715){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  lengthn<-length(sortedx)
  
  if (boot){
    subtract<-t(replicate(times , sort(sample(sortedx, size = 3))))
  }else{
    if (lengthn>300){
      print("Warning: The computational complexity is n^3, bootstrap is recommended.")
    }
    subtract<-t(combn(sortedx, 3))
  }
  getlm<-function(vector){ 
    (1/3)*(vector[3]-2*vector[2]+vector[1])
  }
  dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
  
  getm<-function(vector){ 
    ((1/6)*(2*vector[1]-vector[2]-vector[3])*(-1*vector[1]+2*vector[2]-vector[3])*(-vector[1]-vector[2]+2*vector[3]))
  }
  dp2m<-apply(subtract,MARGIN=1,FUN=getm)
  lm1<-Lmoments(sortedx)
  
  expectdps<-lm1[3]
  expectdp2s<-((sum((sortedx - mean(sortedx))^3)/lengthn)*(lengthn^2/((lengthn-1)*(lengthn-2))))
  dlmo<-mmm(x=dp2lm,interval=9,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
  dmo<-mmm(x=dp2m,interval=9,fast=fast,batch=batch,drm=drm,dqm=dqm)
  
  all<-c(expectdps,dlmo[1],dlmo[2],sd(dp2lm),
         expectdp2s,dmo[1],dmo[2],sd(dp2m))
  return(all)
}
rqtmci<-function(x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.1808,drm=1.7492,dlqm=1.1753,dqm=0.5715,alpha=0.05,nboot=1000){
  data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  pairs1 <- data.frame()
  pairs2 <- data.frame()
  pairs3 <- data.frame()
  pairs4 <- data.frame()
  pairs5 <- data.frame()
  pairs6 <- data.frame()
  for (i in 1:nboot) {
    robustlocation1<-rqtm(data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrm,drm=drm,dlqm=dlqm,dqm=dqm)
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
  estimate=rqtm(x,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrm,drm=drm,dlqm=dlqm,dqm=dqm)
  result <- list(cil3=c(bootlist1[low],bootlist1[up]), cirl3=c(bootlist2[low],bootlist2[up]),ciql3=c(bootlist3[low],bootlist3[up]),
                 citm=c(bootlist4[low],bootlist4[up]),cirtm=c(bootlist5[low],bootlist5[up]),ciqtm=c(bootlist6[low],bootlist6[up]),
                 sel3=sd(bootlist1),serl3=sd(bootlist2),seql3=sd(bootlist3),setm=sd(bootlist4),sertm=sd(bootlist5),seqtm=sd(bootlist6),
                 estimate=estimate)
  return(result)
}

rqfm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=-0.3542,drm=3.4560,dlqm=NaN,dqm=0.1246){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  lengthn<-length(sortedx)
  
  getm<-function(vector){ 
    resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                  4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                  vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                  6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                  vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                               6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
    (resd)
  }
  if (boot){
    subtract<-t(replicate(times , sort(sample(sortedx, size = 4))))
  }else{
    if (lengthn>100){
      print("Warning: The computational complexity is n^4, bootstrap is recommended.")
    }
    subtract<-t(combn(sortedx, 4))
  }
  
  dp2m<-apply(subtract,MARGIN=1,FUN=getm)
  
  getlm<-function(vector){ 
    (1/4)*(vector[4]-3*vector[3]+3*vector[2]-vector[1])
  }

  dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
  
  lm1<-Lmoments(sortedx)
  expectdps<-lm1[4]
  expectdp2s<-(sum((sortedx - mean(sortedx))^4)/lengthn)
  dlmo<-mmm(x=dp2lm,interval=9,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
  dmo<-mmm(x=dp2m,interval=9,fast=fast,batch=batch,drm=drm,dqm=dqm)
  
  all<-c(expectdps,dlmo[1],dlmo[2],sd(dp2lm),
         expectdp2s,dmo[1],dmo[2],sd(dp2m))
  return(all)
}
rqfmci<-function(x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=-0.3542,drm=3.4560,dlqm=NaN,dqm=0.1246,alpha=0.05,nboot=1000){
  data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  pairs1 <- data.frame()
  pairs2 <- data.frame()
  pairs3 <- data.frame()
  pairs4 <- data.frame()
  pairs5 <- data.frame()
  pairs6 <- data.frame()
  for (i in 1:nboot) {
    robustlocation1<-rqfm(data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrm,drm=drm,dlqm=dlqm,dqm=dqm)
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
  estimate=rqfm(x,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrm,drm=drm,dlqm=dlqm,dqm=dqm)
  result <- list(cil4=c(bootlist1[low],bootlist1[up]), cirl4=c(bootlist2[low],bootlist2[up]),ciql4=c(bootlist3[low],bootlist3[up]),
                 cifm=c(bootlist4[low],bootlist4[up]),cirfm=c(bootlist5[low],bootlist5[up]),ciqfm=c(bootlist6[low],bootlist6[up]),
                 sel4=sd(bootlist1),serl4=sd(bootlist2),seql4=sd(bootlist3),sefm=sd(bootlist4),serfm=sd(bootlist5),seqfm=sd(bootlist6),
                 estimate=estimate)
  return(result)
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

NRSssimple<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray")){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  lengthx<-length(sortedx)
  if(standist=="exponential"|| standist=="exp"){
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
    
  }else if (standist=="Rayleigh"|| standist=="Ray"){
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
  mmm1<-mmme(x=sortedx,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  rqscale1<-rqscale(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
  rqtm1<-rqtm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
  rqfm1<-rqfm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
  
  first<-c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4])
  second<-c(l2=rqscale1[1],rl2=rqscale1[2],ql2=rqscale1[3],sdrql2=rqscale1[4],sd=sqrt(rqscale1[5]),
    rsd=sqrt(rqscale1[6]),qsd=sqrt(rqscale1[7]),sdrsd=sqrt(rqscale1[6])*(1/2)*(rqscale1[8]/rqscale1[6]),sdqsd=sqrt(rqscale1[7])*(1/2)*(rqscale1[8]/rqscale1[7])
  )
  third<-c(l3=rqtm1[1]/rqscale1[1],rl3=rqtm1[2]/rqscale1[2],ql3=rqtm1[3]/rqscale1[3],sdrl3=(rqtm1[2]/rqscale1[2])*((rqtm1[4]/rqtm1[2])^2+(rqscale1[4]/rqscale1[2])^2)^(1/2),
            sdql3=(rqtm1[3]/rqscale1[3])*((rqtm1[4]/rqtm1[3])^2+(rqscale1[4]/rqscale1[3])^2)^(1/2),skew=(rqtm1[5])/((rqscale1[5])^(3/2)),rskew=(rqtm1[6])/((rqscale1[6])^(3/2)),qskew=(rqtm1[7])/((rqscale1[7])^(3/2)),
               sdrskew=(rqtm1[6])/((rqscale1[6])^(3/2))*((rqtm1[8]/rqtm1[6])^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2)^(1/2),
               sdqskew=(rqtm1[7])/((rqscale1[7])^(3/2))*((rqtm1[8]/rqtm1[7])^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2)^(1/2)
  )
  fourth<-c(l4=rqfm1[1]/rqscale1[1],rl4=rqfm1[2]/rqscale1[2],ql4=rqfm1[3]/rqscale1[3],sdrl4=(rqfm1[2]/rqscale1[2])*((rqfm1[4]/rqfm1[2])^2+(rqscale1[4]/rqscale1[2])^2)^(1/2),
          sdql4=(rqfm1[3]/rqscale1[3])*((rqfm1[4]/rqfm1[3])^2+(rqscale1[4]/rqscale1[3])^2)^(1/2),kurt=(rqfm1[5])/((rqscale1[5])^(2)),rkurt=(rqfm1[6])/((rqscale1[6])^(2)),qkurt=(rqfm1[7])/((rqscale1[7])^(2)),
    sdrkurt=((rqfm1[6])/((rqscale1[6])^(2)))*((rqfm1[8]/rqfm1[6])^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2)^(1/2),
    sdqkurt=((rqfm1[7])/((rqscale1[7])^(2)))*((rqfm1[8]/rqfm1[7])^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2+((1/2)*(rqscale1[8]/rqscale1[6]))^2+((1/2)*(rqscale1[8]/rqscale1[7]))^2)^(1/2))
  all<-list(first=first,second=second,third=third,fourth=fourth)
  if((rqfm1[7])/((rqscale1[7])^(2))<4.7 & standist=="exponential"){
    print("The quantile kurtosis is lower than 4.7, it might be better to use the Rayleigh distribution as the standard distribution.")
  }
  
  return(all)
}

NRSsci<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray"),alpha=0.05,nboot=100,null_mean=1,null_sd=1,null_skew=2,null_kurt=9,null_l2=0.5,null_l3=1/3,null_l4=1/6){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  lengthx<-length(sortedx)
  if(standist=="exponential"|| standist=="exp"){
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
    
  }else if (standist=="Rayleigh"|| standist=="Ray"){
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
    mmm1<-mmme(x=data[i,],interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
    rqscale1<-rqscale(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
    rqtm1<-rqtm(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
    rqfm1<-rqfm(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
    estimate1<-c(c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4]),
                 
                 c(l2=rqscale1[1],rl2=rqscale1[2],ql2=rqscale1[3]),
                 
                 c(sd=sqrt(rqscale1[5]),rsd=sqrt(rqscale1[6]),qsd=sqrt(rqscale1[7])),
                 
                 c(l3=rqtm1[1]/rqscale1[1],rl3=rqtm1[2]/rqscale1[2],ql3=rqtm1[3]/rqscale1[3]),
                 
                 c(skew=(rqtm1[5])/((rqscale1[5])^(3/2)),rskew=(rqtm1[6])/((rqscale1[6])^(3/2)),qskew=(rqtm1[7])/((rqscale1[7])^(3/2))),
                 
                 c(l4=rqfm1[1]/rqscale1[1],rl4=rqfm1[2]/rqscale1[2],ql4=rqfm1[3]/rqscale1[3]),
                 
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
  pbm<-function(bootlist,null_value){
    p <- mean(bootlist > null_value) + 0.5 * mean(bootlist == null_value)
    p <- 2 * min(c(p, 1 - p))
    p
  }
  mmm1<-mmme(x=sortedx,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  rqscale1<-rqscale(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
  rqtm1<-rqtm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
  rqfm1<-rqfm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
  p_value<-c(mean=pbm(bootlist1,null_mean),etm=pbm(bootlist2,null_mean),rm=pbm(bootlist3,null_mean),qm=pbm(bootlist4,null_mean),
             l2=pbm(bootlist5,null_l2),rl2=pbm(bootlist6,null_l2),ql2=pbm(bootlist7,null_l2),
             sd=pbm(bootlist8,null_sd),
             rsd=pbm(bootlist9,null_sd),qsd=pbm(bootlist10,null_sd),
             l3=pbm(bootlist11,null_l3),rl3=pbm(bootlist12,null_l3),ql3=pbm(bootlist13,null_l3),
             
             skew=pbm(bootlist14,null_skew),rskew=pbm(bootlist15,null_skew),qskew=pbm(bootlist16,null_skew),
             
             l4=pbm(bootlist17,null_l4),rl4=pbm(bootlist18,null_l4),ql4=pbm(bootlist19,null_l4),
             
             kurt=pbm(bootlist20,null_kurt),rkurt=pbm(bootlist21,null_kurt),qkurt=pbm(bootlist22,null_kurt))
  
  
  estimate<-c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4],
              
              l2=rqscale1[1],rl2=rqscale1[2],ql2=rqscale1[3],
              
              sd=sqrt(rqscale1[5]),rsd=sqrt(rqscale1[6]),qsd=sqrt(rqscale1[7]),
              
              l3=rqtm1[1]/rqscale1[1],rl3=rqtm1[2]/rqscale1[2],ql3=rqtm1[3]/rqscale1[3],
              
              skew=(rqtm1[5])/((rqscale1[5])^(3/2)),rskew=(rqtm1[6])/((rqscale1[6])^(3/2)),qskew=(rqtm1[7])/((rqscale1[7])^(3/2)),
              
              l4=rqfm1[1]/rqscale1[1],rl4=rqfm1[2]/rqscale1[2],ql4=rqfm1[3]/rqscale1[3],
              
              kurt=(rqfm1[5])/((rqscale1[5])^(2)),rkurt=(rqfm1[6])/((rqscale1[6])^(2)),qkurt=(rqfm1[7])/((rqscale1[7])^(2)))
  ci<-c(mean=c(bootlist1[low],bootlist1[up]),etm=c(bootlist2[low],bootlist2[up]),rm=c(bootlist3[low],bootlist3[up]),qm=c(bootlist4[low],bootlist4[up]),
        l2=c(bootlist5[low],bootlist5[up]),rl2=c(bootlist6[low],bootlist6[up]),ql2=c(bootlist7[low],bootlist7[up]),
        
        sd=c(bootlist8[low],bootlist8[up]),rsd=c(bootlist9[low],bootlist9[up]),qsd=c(bootlist10[low],bootlist10[up]),
        
        l3=c(bootlist11[low],bootlist11[up]),rl3=c(bootlist12[low],bootlist12[up]),ql3=c(bootlist13[low],bootlist13[up]),
        
        skew=c(bootlist14[low],bootlist14[up]),rskew=c(bootlist15[low],bootlist15[up]),qskew=c(bootlist16[low],bootlist16[up]),
        
        l4=c(bootlist17[low],bootlist17[up]),rl4=c(bootlist18[low],bootlist18[up]),ql4=c(bootlist19[low],bootlist19[up]),
        
        kurt=c(bootlist20[low],bootlist20[up]),rkurt=c(bootlist21[low],bootlist21[up]),qkurt=c(bootlist22[low],bootlist22[up]))
  
  se<-c(mean=sd(bootlist1),etm=sd(bootlist2),rm=sd(bootlist3),qm=sd(bootlist4),
        l2=sd(bootlist5),rl2=sd(bootlist6),ql2=sd(bootlist7),
        
        sd=sd(bootlist8),rsd=sd(bootlist9),qsd=sd(bootlist10),
        
        l3=sd(bootlist11),rl3=sd(bootlist12),ql3=sd(bootlist13),
        
        skew=sd(bootlist14),rskew=sd(bootlist15),qskew=sd(bootlist16),
        
        l4=sd(bootlist17),rl4=sd(bootlist18),ql4=sd(bootlist19),
        
        kurt=sd(bootlist20),rkurt=sd(bootlist21),qkurt=sd(bootlist22))
  
  all<-list(ci=ci,se=se,estimate=estimate,p_value)
  if((rqfm1[7])/((rqscale1[7])^(2))<4.7 & standist=="exponential"){
    print("The quantile kurtosis is lower than 4.7, it might be better to use the Rayleigh distribution as the standard distribution.")
  }
  return(all)
}

NRSsciparallel<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray"),alpha=0.05,nboot=100,null_mean=1,null_sd=1,null_skew=2,null_kurt=9,null_l2=0.5,null_l3=1/3,null_l4=1/6){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  lengthx<-length(sortedx)
  if(standist=="exponential" || standist=="exp"){
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
    
  }else if (standist=="Rayleigh"|| standist=="Ray"){
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
  estimate0<-foreach (i=1:nboot, .combine=rbind) %dopar% {
    library(Lmoments)
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
    
    mmm<-function(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224){
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
      qm1<-quantile(sortedx,quatiletarget)
      rm1<--drm*etm1[3]+etm1[2]+drm*etm1[2]
      names(rm1)<-NULL
      listd<-c(rm1,qm1)
      return(listd)
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
      qm1<-quantile(sortedx,quatiletarget)
      rm1<--drm*etm1[3]+etm1[2]+drm*etm1[2]
      names(rm1)<-NULL
      listd<-c(mean=mean(sortedx),etm=etm1[2],rm=rm1,qm=qm1)
      return(listd)
    }
    
    rqfm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=-0.3542,drm=3.4560,dlqm=NaN,dqm=0.1246){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthn<-length(sortedx)
      
      getm<-function(vector){ 
        resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                      4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                      vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                      6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                      vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                                   6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
        (resd)
      }
      if (boot){
        subtract<-t(replicate(times , sort(sample(sortedx, size = 4))))
      }else{
        if (lengthn>100){
          print("Warning: The computational complexity is n^4, bootstrap is recommended.")
        }
        subtract<-t(combn(sortedx, 4))
      }
      
      dp2m<-apply(subtract,MARGIN=1,FUN=getm)
      
      getlm<-function(vector){ 
        (1/4)*(vector[4]-3*vector[3]+3*vector[2]-vector[1])
      }
      
      dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
      
      lm1<-Lmoments(sortedx)
      expectdps<-lm1[4]
      expectdp2s<-(sum((sortedx - mean(sortedx))^4)/lengthn)
      dlmo<-mmm(x=dp2lm,interval=9,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      dmo<-mmm(x=dp2m,interval=9,fast=fast,batch=batch,drm=drm,dqm=dqm)
      
      all<-c(expectdps,dlmo[1],dlmo[2],sd(dp2lm),
             expectdp2s,dmo[1],dmo[2],sd(dp2m))
      return(all)
    }
    rqscale<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.3659,drm=0.7930,dlqm=0.8218,dqm=0.7825){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      
      lengthn<-length(sortedx)
      if (boot){
        subtract<-t(replicate(times , sort(sample(sortedx, size = 2))))
        getlm<-function(vector){ 
          (vector[2]-vector[1])/2
        }
        dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
        getm<-function(vector){ 
          ((vector[1]-vector[2])^2)/2
        }
        dp2m<-apply(subtract,MARGIN=1,FUN=getm)
      }else{
        if (lengthn>5000){
          print("Warning: The computational complexity is (e*n/2)^2, bootstrap is recommended")
        }
        subtract<-sapply(sortedx, "-", sortedx)
        subtract[lower.tri(subtract)] <- NA
        diag(subtract)=NA
        subtract<-na.omit(as.vector(subtract))
        dp<-subtract[subtract>0]
        dp2lm<-dp/2
        dp2m<-(dp^2)/2
      }
      lm1<-Lmoments(sortedx)
      expectdps<-lm1[2]
      expectdp2s<-(sd(sortedx))^2
      dlmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      dmo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      all<-c(expectdps,dlmo[1],dlmo[2],sd(dp2lm),
             expectdp2s,dmo[1],dmo[2],sd(dp2m)
      )
      return(all)
    }
    rqtm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.1808,drm=1.7492,dlqm=1.1753,dqm=0.5715){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthn<-length(sortedx)
      
      if (boot){
        subtract<-t(replicate(times , sort(sample(sortedx, size = 3))))
      }else{
        if (lengthn>300){
          print("Warning: The computational complexity is n^3, bootstrap is recommended.")
        }
        subtract<-t(combn(sortedx, 3))
      }
      getlm<-function(vector){ 
        (1/3)*(vector[3]-2*vector[2]+vector[1])
      }
      dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
      
      getm<-function(vector){ 
        ((1/6)*(2*vector[1]-vector[2]-vector[3])*(-1*vector[1]+2*vector[2]-vector[3])*(-vector[1]-vector[2]+2*vector[3]))
      }
      dp2m<-apply(subtract,MARGIN=1,FUN=getm)
      lm1<-Lmoments(sortedx)
      
      expectdps<-lm1[3]
      expectdp2s<-((sum((sortedx - mean(sortedx))^3)/lengthn)*(lengthn^2/((lengthn-1)*(lengthn-2))))
      dlmo<-mmm(x=dp2lm,interval=9,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      dmo<-mmm(x=dp2m,interval=9,fast=fast,batch=batch,drm=drm,dqm=dqm)
      
      all<-c(expectdps,dlmo[1],dlmo[2],sd(dp2lm),
             expectdp2s,dmo[1],dmo[2],sd(dp2m))
      return(all)
    }
    mmm1<-mmme(x=data[i,],interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
    rqscale1<-rqscale(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
    rqtm1<-rqtm(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
    rqfm1<-rqfm(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
    estimate1<-c(c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4]),
                 c(l2=rqscale1[1],rl2=rqscale1[2],ql2=rqscale1[3]),
                 c(sd=sqrt(rqscale1[5]),
                   rsd=sqrt(rqscale1[6]),qsd=sqrt(rqscale1[7])),
                 c(l3=rqtm1[1]/rqscale1[1],rl3=rqtm1[2]/rqscale1[2],ql3=rqtm1[3]/rqscale1[3]),
                 
                 c(skew=(rqtm1[5])/((rqscale1[5])^(3/2)),rskew=(rqtm1[6])/((rqscale1[6])^(3/2)),qskew=(rqtm1[7])/((rqscale1[7])^(3/2))),
                 
                 c(l4=rqfm1[1]/rqscale1[1],rl4=rqfm1[2]/rqscale1[2],ql4=rqfm1[3]/rqscale1[3]),
                 
                 c(kurt=(rqfm1[5])/((rqscale1[5])^(2)),rkurt=(rqfm1[6])/((rqscale1[6])^(2)),qkurt=(rqfm1[7])/((rqscale1[7])^(2))))
  }
  
  bootlist1<-sort(as.matrix(estimate0[,1]))
  bootlist2<-sort(as.matrix(estimate0[,2]))
  bootlist3<-sort(as.matrix(estimate0[,3]))
  bootlist4<-sort(as.matrix(estimate0[,4]))
  bootlist5<-sort(as.matrix(estimate0[,5]))
  bootlist6<-sort(as.matrix(estimate0[,6]))
  bootlist7<-sort(as.matrix(estimate0[,7]))
  bootlist8<-sort(as.matrix(estimate0[,8]))
  bootlist9<-sort(as.matrix(estimate0[,9]))
  bootlist10<-sort(as.matrix(estimate0[,10]))
  bootlist11<-sort(as.matrix(estimate0[,11]))
  bootlist12<-sort(as.matrix(estimate0[,12]))
  bootlist13<-sort(as.matrix(estimate0[,13]))
  bootlist14<-sort(as.matrix(estimate0[,14]))
  bootlist15<-sort(as.matrix(estimate0[,15]))
  bootlist16<-sort(as.matrix(estimate0[,16]))
  bootlist17<-sort(as.matrix(estimate0[,17]))
  bootlist18<-sort(as.matrix(estimate0[,18]))
  bootlist19<-sort(as.matrix(estimate0[,19]))
  bootlist20<-sort(as.matrix(estimate0[,20]))
  bootlist21<-sort(as.matrix(estimate0[,21]))
  bootlist22<-sort(as.matrix(estimate0[,22]))
  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  pbm<-function(bootlist,null_value){
    p <- mean(bootlist > null_value) + 0.5 * mean(bootlist == null_value)
    p <- 2 * min(c(p, 1 - p))
    p
  }
  mmm1<-mmme(x=sortedx,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  rqscale1<-rqscale(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
  rqtm1<-rqtm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
  rqfm1<-rqfm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
  p_value<-c(mean=pbm(bootlist1,null_mean),etm=pbm(bootlist2,null_mean),rm=pbm(bootlist3,null_mean),qm=pbm(bootlist4,null_mean),
        l2=pbm(bootlist5,null_l2),rl2=pbm(bootlist6,null_l2),ql2=pbm(bootlist7,null_l2),
        sd=pbm(bootlist8,null_sd),
        rsd=pbm(bootlist9,null_sd),qsd=pbm(bootlist10,null_sd),
        l3=pbm(bootlist11,null_l3),rl3=pbm(bootlist12,null_l3),ql3=pbm(bootlist13,null_l3),
        
        skew=pbm(bootlist14,null_skew),rskew=pbm(bootlist15,null_skew),qskew=pbm(bootlist16,null_skew),
        
        l4=pbm(bootlist17,null_l4),rl4=pbm(bootlist18,null_l4),ql4=pbm(bootlist19,null_l4),
        
        kurt=pbm(bootlist20,null_kurt),rkurt=pbm(bootlist21,null_kurt),qkurt=pbm(bootlist22,null_kurt))
  estimate<-c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4],
              l2=rqscale1[1],rl2=rqscale1[2],ql2=rqscale1[3],
              sd=sqrt(rqscale1[5]),
              rsd=sqrt(rqscale1[6]),qsd=sqrt(rqscale1[7]),
              l3=rqtm1[1]/rqscale1[1],rl3=rqtm1[2]/rqscale1[2],ql3=rqtm1[3]/rqscale1[3],
              
              skew=(rqtm1[5])/((rqscale1[5])^(3/2)),rskew=(rqtm1[6])/((rqscale1[6])^(3/2)),qskew=(rqtm1[7])/((rqscale1[7])^(3/2)),
              
              l4=rqfm1[1]/rqscale1[1],rl4=rqfm1[2]/rqscale1[2],ql4=rqfm1[3]/rqscale1[3],
              
              kurt=(rqfm1[5])/((rqscale1[5])^(2)),rkurt=(rqfm1[6])/((rqscale1[6])^(2)),qkurt=(rqfm1[7])/((rqscale1[7])^(2)))
  ci<-c(mean=c(bootlist1[low],bootlist1[up]),etm=c(bootlist2[low],bootlist2[up]),rm=c(bootlist3[low],bootlist3[up]),qm=c(bootlist4[low],bootlist4[up]),
        l2=c(bootlist5[low],bootlist5[up]),rl2=c(bootlist6[low],bootlist6[up]),ql2=c(bootlist7[low],bootlist7[up]),
        sd=c(bootlist8[low],bootlist8[up]),
        rsd=c(bootlist9[low],bootlist9[up]),qsd=c(bootlist10[low],bootlist10[up]),
        l3=c(bootlist11[low],bootlist11[up]),rl3=c(bootlist12[low],bootlist12[up]),ql3=c(bootlist13[low],bootlist13[up]),
        
        skew=c(bootlist14[low],bootlist14[up]),rskew=c(bootlist15[low],bootlist15[up]),qskew=c(bootlist16[low],bootlist16[up]),
        
        l4=c(bootlist17[low],bootlist17[up]),rl4=c(bootlist18[low],bootlist18[up]),ql4=c(bootlist19[low],bootlist19[up]),
        
        kurt=c(bootlist20[low],bootlist20[up]),rkurt=c(bootlist21[low],bootlist21[up]),qkurt=c(bootlist22[low],bootlist22[up]))
  
  se<-c(mean=sd(bootlist1),etm=sd(bootlist2),rm=sd(bootlist3),qm=sd(bootlist4),
        l2=sd(bootlist5),rl2=sd(bootlist6),ql2=sd(bootlist7),
        sd=sd(bootlist8),
        rsd=sd(bootlist9),qsd=sd(bootlist10),
        l3=sd(bootlist11),rl3=sd(bootlist12),ql3=sd(bootlist13),
        
        skew=sd(bootlist14),rskew=sd(bootlist15),qskew=sd(bootlist16),
        
        l4=sd(bootlist17),rl4=sd(bootlist18),ql4=sd(bootlist19),
        
        kurt=sd(bootlist20),rkurt=sd(bootlist21),qkurt=sd(bootlist22))
  
  all<-list(ci=ci,se=se,estimate=estimate,p_value=p_value)
  if((rqfm1[7])/((rqscale1[7])^(2))<4.7 & standist=="exponential"){
    print("The quantile kurtosis is lower than 4.7, it might be better to use the Rayleigh distribution as the standard distribution.")
  }
  return(all)
}
pbh2parallel<-function (x,y,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray"),alpha=0.05,nboot=100){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  lengthx<-length(sortedx)
  sortedy<-sort(y,decreasing = FALSE,method ="radix")
  lengthy<-length(sortedy)
  if(standist=="exponential" || standist=="exp"){
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
    
  }else if (standist=="Rayleigh"|| standist=="Ray"){
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
  datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
  
  estimate0<-foreach (i=1:nboot, .combine=rbind) %dopar% {
    library(Lmoments)
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
    
    mmm<-function(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224){
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
      qm1<-quantile(sortedx,quatiletarget)
      rm1<--drm*etm1[3]+etm1[2]+drm*etm1[2]
      names(rm1)<-NULL
      listd<-c(rm1,qm1)
      return(listd)
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
      qm1<-quantile(sortedx,quatiletarget)
      rm1<--drm*etm1[3]+etm1[2]+drm*etm1[2]
      names(rm1)<-NULL
      listd<-c(mean=mean(sortedx),etm=etm1[2],rm=rm1,qm=qm1)
      return(listd)
    }
    
    rqfm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=-0.3542,drm=3.4560,dlqm=NaN,dqm=0.1246){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthn<-length(sortedx)
      
      getm<-function(vector){ 
        resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                      4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                      vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                      6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                      vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                                   6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
        (resd)
      }
      if (boot){
        subtract<-t(replicate(times , sort(sample(sortedx, size = 4))))
      }else{
        if (lengthn>100){
          print("Warning: The computational complexity is n^4, bootstrap is recommended.")
        }
        subtract<-t(combn(sortedx, 4))
      }
      
      dp2m<-apply(subtract,MARGIN=1,FUN=getm)
      
      getlm<-function(vector){ 
        (1/4)*(vector[4]-3*vector[3]+3*vector[2]-vector[1])
      }
      
      dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
      
      lm1<-Lmoments(sortedx)
      expectdps<-lm1[4]
      expectdp2s<-(sum((sortedx - mean(sortedx))^4)/lengthn)
      dlmo<-mmm(x=dp2lm,interval=9,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      dmo<-mmm(x=dp2m,interval=9,fast=fast,batch=batch,drm=drm,dqm=dqm)
      
      all<-c(expectdps,dlmo[1],dlmo[2],sd(dp2lm),
             expectdp2s,dmo[1],dmo[2],sd(dp2m))
      return(all)
    }
    rqscale<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.3659,drm=0.7930,dlqm=0.8218,dqm=0.7825){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      
      lengthn<-length(sortedx)
      if (boot){
        subtract<-t(replicate(times , sort(sample(sortedx, size = 2))))
        getlm<-function(vector){ 
          (vector[2]-vector[1])/2
        }
        dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
        getm<-function(vector){ 
          ((vector[1]-vector[2])^2)/2
        }
        dp2m<-apply(subtract,MARGIN=1,FUN=getm)
      }else{
        if (lengthn>5000){
          print("Warning: The computational complexity is (e*n/2)^2, bootstrap is recommended")
        }
        subtract<-sapply(sortedx, "-", sortedx)
        subtract[lower.tri(subtract)] <- NA
        diag(subtract)=NA
        subtract<-na.omit(as.vector(subtract))
        dp<-subtract[subtract>0]
        dp2lm<-dp/2
        dp2m<-(dp^2)/2
      }
      lm1<-Lmoments(sortedx)
      expectdps<-lm1[2]
      expectdp2s<-(sd(sortedx))^2
      dlmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      dmo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      all<-c(expectdps,dlmo[1],dlmo[2],sd(dp2lm),
             expectdp2s,dmo[1],dmo[2],sd(dp2m)
      )
      return(all)
    }
    rqtm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.1808,drm=1.7492,dlqm=1.1753,dqm=0.5715){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthn<-length(sortedx)
      
      if (boot){
        subtract<-t(replicate(times , sort(sample(sortedx, size = 3))))
      }else{
        if (lengthn>300){
          print("Warning: The computational complexity is n^3, bootstrap is recommended.")
        }
        subtract<-t(combn(sortedx, 3))
      }
      getlm<-function(vector){ 
        (1/3)*(vector[3]-2*vector[2]+vector[1])
      }
      dp2lm<-apply(subtract,MARGIN=1,FUN=getlm)
      
      getm<-function(vector){ 
        ((1/6)*(2*vector[1]-vector[2]-vector[3])*(-1*vector[1]+2*vector[2]-vector[3])*(-vector[1]-vector[2]+2*vector[3]))
      }
      dp2m<-apply(subtract,MARGIN=1,FUN=getm)
      lm1<-Lmoments(sortedx)
      
      expectdps<-lm1[3]
      expectdp2s<-((sum((sortedx - mean(sortedx))^3)/lengthn)*(lengthn^2/((lengthn-1)*(lengthn-2))))
      dlmo<-mmm(x=dp2lm,interval=9,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      dmo<-mmm(x=dp2m,interval=9,fast=fast,batch=batch,drm=drm,dqm=dqm)
      
      all<-c(expectdps,dlmo[1],dlmo[2],sd(dp2lm),
             expectdp2s,dmo[1],dmo[2],sd(dp2m))
      return(all)
    }
    mmmx<-mmme(x=datax[i,],interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
    rqscalex<-rqscale(x=datax[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
    rqtmx<-rqtm(x=datax[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
    rqfmx<-rqfm(x=datax[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
    mmmy<-mmme(x=datay[i,],interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
    rqscaley<-rqscale(x=datay[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
    rqtmy<-rqtm(x=datay[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
    rqfmy<-rqfm(x=datay[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
    estimate1<-c(c(meanx=mmmx[1],etmx=mmmx[2],rmx=mmmx[3],qmx=mmmx[4]),
                 c(l2x=rqscalex[1],rl2x=rqscalex[2],ql2x=rqscalex[3]),
                 c(sdx=sqrt(rqscalex[5]),
                   rsdx=sqrt(rqscalex[6]),qsdx=sqrt(rqscalex[7])),
                 c(l3x=rqtmx[1]/rqscalex[1],rl3x=rqtmx[2]/rqscalex[2],ql3x=rqtmx[3]/rqscalex[3]),
                 
                 c(skewx=(rqtmx[5])/((rqscalex[5])^(3/2)),rskewx=(rqtmx[6])/((rqscalex[6])^(3/2)),qskewx=(rqtmx[7])/((rqscalex[7])^(3/2))),
                 
                 c(l4x=rqfmx[1]/rqscalex[1],rl4x=rqfmx[2]/rqscalex[2],ql4x=rqfmx[3]/rqscalex[3]),
                 
                 c(kurtx=(rqfmx[5])/((rqscalex[5])^(2)),rkurtx=(rqfmx[6])/((rqscalex[6])^(2)),qkurtx=(rqfmx[7])/((rqscalex[7])^(2))),
                 
                 c(meany=mmmy[1],etmy=mmmy[2],rmy=mmmy[3],qmy=mmmy[4]),
                 c(l2y=rqscaley[1],rl2y=rqscaley[2],ql2y=rqscaley[3]),
                 c(sdy=sqrt(rqscaley[5]),
                   rsdy=sqrt(rqscaley[6]),qsd=sqrt(rqscaley[7])),
                 c(l3y=rqtmy[1]/rqscaley[1],rl3y=rqtmy[2]/rqscaley[2],ql3y=rqtmy[3]/rqscaley[3]),
                 
                 c(skewy=(rqtmy[5])/((rqscaley[5])^(3/2)),rskewy=(rqtmy[6])/((rqscaley[6])^(3/2)),qskewy=(rqtmy[7])/((rqscaley[7])^(3/2))),
                 
                 c(l4y=rqfmy[1]/rqscaley[1],rl4y=rqfmy[2]/rqscaley[2],ql4y=rqfmy[3]/rqscaley[3]),
                 
                 c(kurty=(rqfmy[5])/((rqscaley[5])^(2)),rkurt=(rqfmy[6])/((rqscaley[6])^(2)),qkurty=(rqfmy[7])/((rqscaley[7])^(2)))
                 )
  }
  
  bootlist1<-(as.matrix(estimate0[,1]))
  bootlist2<-(as.matrix(estimate0[,2]))
  bootlist3<-(as.matrix(estimate0[,3]))
  bootlist4<-(as.matrix(estimate0[,4]))
  bootlist5<-(as.matrix(estimate0[,5]))
  bootlist6<-(as.matrix(estimate0[,6]))
  bootlist7<-(as.matrix(estimate0[,7]))
  bootlist8<-(as.matrix(estimate0[,8]))
  bootlist9<-(as.matrix(estimate0[,9]))
  bootlist10<-(as.matrix(estimate0[,10]))
  bootlist11<-(as.matrix(estimate0[,11]))
  bootlist12<-(as.matrix(estimate0[,12]))
  bootlist13<-(as.matrix(estimate0[,13]))
  bootlist14<-(as.matrix(estimate0[,14]))
  bootlist15<-(as.matrix(estimate0[,15]))
  bootlist16<-(as.matrix(estimate0[,16]))
  bootlist17<-(as.matrix(estimate0[,17]))
  bootlist18<-(as.matrix(estimate0[,18]))
  bootlist19<-(as.matrix(estimate0[,19]))
  bootlist20<-(as.matrix(estimate0[,20]))
  bootlist21<-(as.matrix(estimate0[,21]))
  bootlist22<-(as.matrix(estimate0[,22]))
  
  bootlist23<-(as.matrix(estimate0[,23]))
  bootlist24<-(as.matrix(estimate0[,24]))
  bootlist25<-(as.matrix(estimate0[,25]))
  bootlist26<-(as.matrix(estimate0[,26]))
  bootlist27<-(as.matrix(estimate0[,27]))
  bootlist28<-(as.matrix(estimate0[,28]))
  bootlist29<-(as.matrix(estimate0[,29]))
  bootlist30<-(as.matrix(estimate0[,30]))
  bootlist31<-(as.matrix(estimate0[,31]))
  bootlist32<-(as.matrix(estimate0[,32]))
  bootlist33<-(as.matrix(estimate0[,33]))
  bootlist34<-(as.matrix(estimate0[,34]))
  bootlist35<-(as.matrix(estimate0[,35]))
  bootlist36<-(as.matrix(estimate0[,36]))
  bootlist37<-(as.matrix(estimate0[,37]))
  bootlist38<-(as.matrix(estimate0[,38]))
  bootlist39<-(as.matrix(estimate0[,39]))
  bootlist40<-(as.matrix(estimate0[,40]))
  bootlist41<-(as.matrix(estimate0[,41]))
  bootlist42<-(as.matrix(estimate0[,42]))
  bootlist43<-(as.matrix(estimate0[,43]))
  bootlist44<-(as.matrix(estimate0[,44]))
  
  bootlist110<-sort(bootlist1-bootlist23)
  bootlist220<-sort(bootlist2-bootlist24)
  bootlist330<-sort(bootlist3-bootlist25)
  bootlist440<-sort(bootlist4-bootlist26)
  bootlist550<-sort(bootlist5-bootlist27)
  bootlist660<-sort(bootlist6-bootlist28)
  bootlist770<-sort(bootlist7-bootlist29)
  bootlist880<-sort(bootlist8-bootlist30)
  bootlist990<-sort(bootlist9-bootlist31)
  bootlist1010<-sort(bootlist10-bootlist32)
  bootlist1111<-sort(bootlist11-bootlist33)
  bootlist1212<-sort(bootlist12-bootlist34)
  bootlist1313<-sort(bootlist13-bootlist35)
  bootlist1414<-sort(bootlist14-bootlist36)
  bootlist1515<-sort(bootlist15-bootlist37)
  bootlist1616<-sort(bootlist16-bootlist38)
  bootlist1717<-sort(bootlist17-bootlist39)
  bootlist1818<-sort(bootlist18-bootlist40)
  bootlist1919<-sort(bootlist19-bootlist41)
  bootlist2020<-sort(bootlist20-bootlist42)
  bootlist2121<-sort(bootlist21-bootlist43)
  bootlist2222<-sort(bootlist22-bootlist44)
  
  
  
  bootlist1<-sort(bootlist1)
  bootlist2<-sort(bootlist2)
  bootlist3<-sort(bootlist3)
  bootlist4<-sort(bootlist4)
  bootlist5<-sort(bootlist5)
  bootlist6<-sort(bootlist6)
  bootlist7<-sort(bootlist7)
  bootlist8<-sort(bootlist8)
  bootlist9<-sort(bootlist9)
  bootlist10<-sort(bootlist10)
  bootlist11<-sort(bootlist11)
  bootlist12<-sort(bootlist12)
  bootlist13<-sort(bootlist13)
  bootlist14<-sort(bootlist14)
  bootlist15<-sort(bootlist15)
  bootlist16<-sort(bootlist16)
  bootlist17<-sort(bootlist17)
  bootlist18<-sort(bootlist18)
  bootlist19<-sort(bootlist19)
  bootlist20<-sort(bootlist20)
  bootlist21<-sort(bootlist21)
  bootlist22<-sort(bootlist22)
  
  bootlist23<-sort(bootlist23)
  bootlist24<-sort(bootlist24)
  bootlist25<-sort(bootlist25)
  bootlist26<-sort(bootlist26)
  bootlist27<-sort(bootlist27)
  bootlist28<-sort(bootlist28)
  bootlist29<-sort(bootlist29)
  bootlist30<-sort(bootlist30)
  bootlist31<-sort(bootlist31)
  bootlist32<-sort(bootlist32)
  bootlist33<-sort(bootlist33)
  bootlist34<-sort(bootlist34)
  bootlist35<-sort(bootlist35)
  bootlist36<-sort(bootlist36)
  bootlist37<-sort(bootlist37)
  bootlist38<-sort(bootlist38)
  bootlist39<-sort(bootlist39)
  bootlist40<-sort(bootlist40)
  bootlist41<-sort(bootlist41)
  bootlist42<-sort(bootlist42)
  bootlist43<-sort(bootlist43)
  bootlist44<-sort(bootlist44)
  
  
  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  
  pb2<-function(bootlist){
    p <- mean(bootlist < 0) + 0.5 * mean(bootlist == 0)
    p <- 2 * min(c(p, 1 - p))
    p
  }
  ci_diff<-c(mean=c(bootlist110[low],bootlist110[up]),etm=c(bootlist220[low],bootlist220[up]),rm=c(bootlist330[low],bootlist330[up]),qm=c(bootlist440[low],bootlist440[up]),
        l2=c(bootlist550[low],bootlist550[up]),rl2=c(bootlist660[low],bootlist660[up]),ql2=c(bootlist770[low],bootlist770[up]),
        sd=c(bootlist880[low],bootlist880[up]),
        rsd=c(bootlist990[low],bootlist990[up]),qsd=c(bootlist1010[low],bootlist1010[up]),
        l3=c(bootlist1111[low],bootlist1111[up]),rl3=c(bootlist1212[low],bootlist1212[up]),ql3=c(bootlist1313[low],bootlist1313[up]),
        
        skew=c(bootlist1414[low],bootlist1414[up]),rskew=c(bootlist1515[low],bootlist1515[up]),qskew=c(bootlist1616[low],bootlist1616[up]),
        
        l4=c(bootlist1717[low],bootlist1717[up]),rl4=c(bootlist1818[low],bootlist1818[up]),ql4=c(bootlist1919[low],bootlist1919[up]),
        
        kurt=c(bootlist2020[low],bootlist2020[up]),rkurt=c(bootlist2121[low],bootlist2121[up]),qkurt=c(bootlist2222[low],bootlist2222[up]))
  
  p_value_diff<-c(mean=pb2(bootlist110),etm=pb2(bootlist220),rm=pb2(bootlist330),qm=pb2(bootlist440),
             l2=pb2(bootlist550),rl2=pb2(bootlist660),ql2=pb2(bootlist770),
             sd=pb2(bootlist880),
             rsd=pb2(bootlist990),qsd=pb2(bootlist1010),
             l3=pb2(bootlist1111),rl3=pb2(bootlist1212),ql3=pb2(bootlist1313),
             
             skew=pb2(bootlist1414),rskew=pb2(bootlist1515),qskew=pb2(bootlist1616),
             
             l4=pb2(bootlist1717),rl4=pb2(bootlist1818),ql4=pb2(bootlist1919),
             
             kurt=pb2(bootlist2020),rkurt=pb2(bootlist2121),qkurt=pb2(bootlist2222))
  
  mmmx<-mmme(x=sortedx,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  rqscalex<-rqscale(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
  rqtmx<-rqtm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
  rqfmx<-rqfm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
  mmmy<-mmme(x=sortedy,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  rqscaley<-rqscale(x=sortedy,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale)
  rqtmy<-rqtm(x=sortedy,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm)
  rqfmy<-rqfm(x=sortedy,interval=interval,fast=fast,batch=batch,boot=boot,times =times ,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm)
  estimate<-c(c(meanx=mmmx[1],etmx=mmmx[2],rmx=mmmx[3],qmx=mmmx[4]),
               c(l2x=rqscalex[1],rl2x=rqscalex[2],ql2x=rqscalex[3]),
               c(sdx=sqrt(rqscalex[5]),
                 rsdx=sqrt(rqscalex[6]),qsdx=sqrt(rqscalex[7])),
               c(l3x=rqtmx[1]/rqscalex[1],rl3x=rqtmx[2]/rqscalex[2],ql3x=rqtmx[3]/rqscalex[3]),
               
               c(skewx=(rqtmx[5])/((rqscalex[5])^(3/2)),rskewx=(rqtmx[6])/((rqscalex[6])^(3/2)),qskewx=(rqtmx[7])/((rqscalex[7])^(3/2))),
               
               c(l4x=rqfmx[1]/rqscalex[1],rl4x=rqfmx[2]/rqscalex[2],ql4x=rqfmx[3]/rqscalex[3]),
               
               c(kurtx=(rqfmx[5])/((rqscalex[5])^(2)),rkurtx=(rqfmx[6])/((rqscalex[6])^(2)),qkurtx=(rqfmx[7])/((rqscalex[7])^(2))),
               
               c(meany=mmmy[1],etmy=mmmy[2],rmy=mmmy[3],qmy=mmmy[4]),
               c(l2y=rqscaley[1],rl2y=rqscaley[2],ql2y=rqscaley[3]),
               c(sdy=sqrt(rqscaley[5]),
                 rsdy=sqrt(rqscaley[6]),qsd=sqrt(rqscaley[7])),
               c(l3y=rqtmy[1]/rqscaley[1],rl3y=rqtmy[2]/rqscaley[2],ql3y=rqtmy[3]/rqscaley[3]),
               
               c(skewy=(rqtmy[5])/((rqscaley[5])^(3/2)),rskewy=(rqtmy[6])/((rqscaley[6])^(3/2)),qskewy=(rqtmy[7])/((rqscaley[7])^(3/2))),
               
               c(l4y=rqfmy[1]/rqscaley[1],rl4y=rqfmy[2]/rqscaley[2],ql4y=rqfmy[3]/rqscaley[3]),
               
               c(kurty=(rqfmy[5])/((rqscaley[5])^(2)),rkurt=(rqfmy[6])/((rqscaley[6])^(2)),qkurty=(rqfmy[7])/((rqscaley[7])^(2)))
  )
  
  ci<-c(meanx=c(bootlist1[low],bootlist1[up]),etmx=c(bootlist2[low],bootlist2[up]),rmx=c(bootlist3[low],bootlist3[up]),qmx=c(bootlist4[low],bootlist4[up]),
        l2x=c(bootlist5[low],bootlist5[up]),rl2x=c(bootlist6[low],bootlist6[up]),ql2x=c(bootlist7[low],bootlist7[up]),
        sdx=c(bootlist8[low],bootlist8[up]),
        rsdx=c(bootlist9[low],bootlist9[up]),qsdx=c(bootlist10[low],bootlist10[up]),
        l3x=c(bootlist11[low],bootlist11[up]),rl3x=c(bootlist12[low],bootlist12[up]),ql3x=c(bootlist13[low],bootlist13[up]),
        
        skewx=c(bootlist14[low],bootlist14[up]),rskewx=c(bootlist15[low],bootlist15[up]),qskewx=c(bootlist16[low],bootlist16[up]),
        
        l4x=c(bootlist17[low],bootlist17[up]),rl4x=c(bootlist18[low],bootlist18[up]),ql4x=c(bootlist19[low],bootlist19[up]),
        
        kurtx=c(bootlist20[low],bootlist20[up]),rkurtx=c(bootlist21[low],bootlist21[up]),qkurtx=c(bootlist22[low],bootlist22[up]),
        
        meany=c(bootlist23[low],bootlist23[up]),etmy=c(bootlist24[low],bootlist24[up]),rmy=c(bootlist25[low],bootlist25[up]),qmy=c(bootlist26[low],bootlist26[up]),
        l2y=c(bootlist27[low],bootlist27[up]),rl2y=c(bootlist28[low],bootlist28[up]),ql2y=c(bootlist29[low],bootlist29[up]),
        sdy=c(bootlist30[low],bootlist30[up]),
        rsdy=c(bootlist31[low],bootlist31[up]),qsdy=c(bootlist32[low],bootlist32[up]),
        l3y=c(bootlist33[low],bootlist33[up]),rl3y=c(bootlist34[low],bootlist34[up]),ql3y=c(bootlist35[low],bootlist35[up]),
        
        skewy=c(bootlist36[low],bootlist36[up]),rskewy=c(bootlist37[low],bootlist37[up]),qskewy=c(bootlist38[low],bootlist38[up]),
        
        l4y=c(bootlist39[low],bootlist39[up]),rl4y=c(bootlist40[low],bootlist40[up]),ql4y=c(bootlist41[low],bootlist41[up]),
        
        kurty=c(bootlist42[low],bootlist42[up]),rkurty=c(bootlist43[low],bootlist43[up]),qkurty=c(bootlist44[low],bootlist44[up])
        )
  

  se<-c(meanx=sd(bootlist1),etmx=sd(bootlist2),rmx=sd(bootlist3),qmx=sd(bootlist4),
        l2x=sd(bootlist5),rl2x=sd(bootlist6),ql2x=sd(bootlist7),
        sdx=sd(bootlist8),
        rsdx=sd(bootlist9),qsdx=sd(bootlist10),
        l3x=sd(bootlist11),rl3x=sd(bootlist12),ql3x=sd(bootlist13),
        
        skewx=sd(bootlist14),rskewx=sd(bootlist15),qskewx=sd(bootlist16),
        
        l4x=sd(bootlist17),rl4x=sd(bootlist18),ql4x=sd(bootlist19),
        
        kurtx=sd(bootlist20),rkurtx=sd(bootlist21),qkurtx=sd(bootlist22),
        
        meany=sd(bootlist23),etmy=sd(bootlist24),rmy=sd(bootlist25),qmy=sd(bootlist26),
        l2y=sd(bootlist27),rl2y=sd(bootlist28),ql2y=sd(bootlist29),
        sdy=sd(bootlist30),
        rsdy=sd(bootlist31),qsdy=sd(bootlist32),
        l3y=sd(bootlist33),rl3y=sd(bootlist34),ql3y=sd(bootlist35),
        
        skewy=sd(bootlist36),rskewy=sd(bootlist37),qskewy=sd(bootlist38),
        
        l4y=sd(bootlist39),rl4y=sd(bootlist40),ql4y=sd(bootlist41),
        
        kurty=sd(bootlist42),rkurty=sd(bootlist43),qkurty=sd(bootlist44)
        )
  
  all<-list(p_value_diff,ci_diff,ci=ci,se=se,estimate=estimate)
  if((rqfmx[7])/((rqscalex[7])^(2))<4.7 & standist=="exponential"){
    print("The quantile kurtosis is lower than 4.7, it might be better to use the Rayleigh distribution as the standard distribution.")
  }
  return(all)
}


NRSs<-function(x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray"),cise = FALSE,parallel=TRUE,alpha = 0.05,nboot = 100,null_mean=1,null_sd=1,null_skew=2,null_kurt=9,null_l2=0.5,null_l3=1/3,null_l4=1/6){
  if (times%%9!=0){
    return ("Please set times as a multiple of 9.")
  }
  if(cise & parallel){
    return (NRSsciparallel(x, interval=interval,fast=fast,batch=batch,boot=boot,times =times ,standist=standist,alpha=alpha,nboot=nboot,null_mean=null_mean,null_sd=null_sd,null_skew=null_skew,null_kurt=null_kurt,null_l2=null_l2,null_l3=null_l3,null_l4=null_l4))
  } else if(cise){return (NRSsci(x, interval=interval,fast=fast,batch=batch,boot=boot,times =times ,standist=standist,alpha=alpha,nboot=nboot,null_mean=null_mean,null_sd=null_sd,null_skew=null_skew,null_kurt=null_kurt,null_l2=null_l2,null_l3=null_l3,null_l4=null_l4))
  }
  else{return (NRSssimple(x, interval=interval,fast=fast,batch=batch,boot=boot,times =times ,standist=standist))
}}

#test
xexp<-rexp(5400,1)

#the population mean is 1
#the population standard deviation is 1
#the population L2-moment is 1/2
#the population skewness is 2
#the population L-skewness is 1/3
#the population kurtosis is 9
#the population L-kurtosis is 1/6
#no d for ql4, because the distribution of U-statistic of L4-moment does not follow mean-ETM-median inequality

#this standard deviation of the distribution of U-statistic is calculated based on the law of prorogation of uncertainty.
NRSs(x=xexp,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exp",cise = FALSE,parallel=TRUE,alpha = 0.05,nboot = 100)
#Arguments
#x:a numeric vector
#interval: The b value in equinterval trimmed mean and complement trimmed mean, notifying that the breakdown points for higher order moments/L-moments are b*k, not b.
#fast: logical; if "TRUE", the approximation solution based on data augmentation for n mod b =/ 0 is used, only available when the sample size is smaller than 10000.
#batch: if fast=FALSE, the approximation solution based on multiple-imputation, the "auto" option is 500000/(length of x), which is corresponding to five decimal accuracy.
#boot: logical; if "TRUE", bootstrap is used for second and higher order moments/L-moments estimations, if not used, the computational time is often unacceptable due to the combinatorial explosion.
#times : the number of subsampling times, a multiple of 9, used in bootstrap.
#standist: a character string giving the standard distribution to be used to calibrate the d value. This must partially match either "exponential" or "Rayleigh", with default "exponential" and may be abbreviated to a unique prefix (the first three letters).
#cise: logical; if "TRUE", the confidence interval and standard error will be estimated using bootstrap.  
#parallel: logical; whether use parallel computing for the confidential interval and standard error, if not used, 100 nboot takes >10 mins, while if used, in a typical PC, the running time is about 1 min. Additional foreach and doparallel packages are required.
#alpha: the alpha level for confidence interval computation.
#nboot: the number of bootstrap samples for confidence interval computation.

#To make comparisons easier, sample standardized moments and scaled L-moments are provided. 
#The standard deviations of the distributions of U-statistic can be used to estimate the consistency percentage.
#The standard error and confidential interval of the robust or quantile mean can be accurately estimated by bootstrapping.

rqmean(x=xexp,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = TRUE,alpha = 0.05,nboot = 1000)
#A similar approach can be applied to all NRSs, but just 100 nboot takes ~10 mins.

#If you don't want to wait for a long time to compare, don't run the following code.
#NRSs(x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,parallel=FALSE,alpha = 0.05,nboot = 100)

#A solution is parallel computing (takes 1 min with 12 cores, but is unavailable on some types of computers).

#The standard errors of robust skewness and kurtosis are lower than those of sample skewness and kurtosis.

#also, the one-sample hypothesis testing can be done with percentile bootstrap method.

#the null values are the corresponding population parameters

NRSs(x=xexp,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,parallel=TRUE,alpha = 0.05,nboot=100,null_mean=1,null_sd=1,null_skew=2,null_kurt=9,null_l2=0.5,null_l3=1/3,null_l4=1/6)



#two-goup comparison can also be done with a similar approach.
xexp<-rexp(5400,1.1)
yexp<-rexp(5400,1)

#to reduce the test time, the boot times of U-statistics are 5400, instead of 54000. if 54000, takes around 3 mins.
pbh2parallel(x=xexp,y=yexp,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

xexp<-c(rexp(5380,1),rnorm(20,10))
yexp<-rexp(5400,1)

#test of outliers
pbh2parallel(x=xexp,y=yexp,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

xexp<-c(rexp(5400,1))
yexp<-rexp(5400,1)

#test of null hypothesis
pbh2parallel(x=xexp,y=yexp,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)


#It should be noticed that, the percentile bootstrap method is appealing. John Rice, Mathematical Statistics and Data Analysis, 2nd edition, p. 272
#so this method should be further testing....

library(lmom)


a=500
xgamma<-c(rgamma(5400, shape=a/100, rate = 1))
targetgam<-lmrgam(para = c(a/100, 1), nmom = 4)
NRSs(x=xgamma,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
#A rule of thumb for desired consistency performance (all four moments > 90%) is that 
#the kurtosis of the underlying distribution should be within [1/2,2] times that of the standard distribution used to calibrate the d values. 
#That means, using exponential as the standard distribution, the kurtosis should be within 4.5 to 18.

#While accurately estimating population kurtosis is hard, finding a rough range and choosing the right standard should be easy in practice.
#Also, if the quantile kurtosis is less than 4.7, it highly indicates the need to change to Rayleigh.
NRSs(x=xgamma,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
#parallel is default for confidential interval.
NRSs(x=xgamma,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=targetgam[1],null_sd=sqrt((a/100)),null_skew=2/sqrt(a/100),null_kurt=((6/(a/100))+3),null_l2=targetgam[2],null_l3=targetgam[3],null_l4=targetgam[4])

NRSs(x=xgamma,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=targetgam[1],null_sd=sqrt((a/100)),null_skew=2/sqrt(a/100),null_kurt=((6/(a/100))+3),null_l2=targetgam[2],null_l3=targetgam[3],null_l4=targetgam[4])

a=150
xgamma<-c(rgamma(5400, shape=a/100, rate = 1))
targetgam<-lmrgam(para = c(a/100, 1), nmom = 4)
NRSs(x=xgamma,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xgamma,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
#even the kurtosis is not very high, 7, the standard errors are still lower than sample moments. 
NRSs(x=xgamma,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=targetgam[1],null_sd=sqrt((a/100)),null_skew=2/sqrt(a/100),null_kurt=((6/(a/100))+3),null_l2=targetgam[2],null_l3=targetgam[3],null_l4=targetgam[4])




xRayleigh<-rRayleigh(n=5400, scale = 1) 

NRSs(x=xRayleigh,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xRayleigh,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
NRSs(x=xRayleigh,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh",cise = TRUE,parallel=TRUE,alpha = 0.05,nboot=100,null_mean=sqrt(pi/2),null_sd=sqrt(2-(pi/2)),null_skew=2*sqrt((pi))*(pi-3)/((4-pi)^(3/2)),null_kurt=(3-(6*(pi)^2-24*(pi)+16)/((4-pi)^(2))),null_l2=0.5*(sqrt(2)-1)*sqrt(pi),null_l3=((1/6)*(2*sqrt(6)+3*sqrt(2)-9)*sqrt(pi))/(0.5*(sqrt(2)-1)*sqrt(pi)),null_l4=((sqrt((77/6)-5*sqrt(6))-3/4)*sqrt(2*pi))/(0.5*(sqrt(2)-1)*sqrt(pi)))


xRayleigh<-rRayleigh(n=5400, scale = 1.08) 
yRayleigh<-rRayleigh(n=5400, scale = 1) 
pbh2parallel(x=xRayleigh,y=yRayleigh,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)




xpois<-rpois(5400,8)
#the population mean is 8
#quantile mean returns 7 or 9
#the population variance is 8
#the population L2-moment is 1.583
#the population skewness is 0.3535534
#the population L-skewness is 0.0592
#the population kurtosis is 1/8+3
#the population L-kurtosis is 0.1204

#quantile L-moments are not suitable for discrete distributions

#because the kurtosis of poisson is close to 3, the biases of robust/quantile kurtosis based on the exponential distribution are large.

NRSs(x=xpois,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")

NRSs(x=xpois,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
NRSs(x=xpois,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh",cise = TRUE,parallel=TRUE,alpha = 0.05,nboot=100,null_mean=8,null_sd=sqrt(8),null_skew=0.3535534,null_kurt=(1/8+3),null_l2=1.583,null_l3=0.0592,null_l4=0.1204)


xnorm<-c(rnorm(5400))

NRSs(x=xnorm,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xnorm,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")

NRSs(x=xnorm,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=0,null_sd=1,null_skew=0,null_kurt=3,null_l2=1/sqrt(pi),null_l3=0,null_l4=(30*(1/(pi))*(atan(sqrt(2)))-9))

xlogis<-c(rlogis(5400, location = 0, scale = 1))

NRSs(x=xlogis,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xlogis,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
NRSs(x=xlogis,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=0,null_sd=sqrt(((pi^2)/3)),null_skew=0,null_kurt=(((6/5)+3)*((sqrt((pi^2)/3))^4))/(sqrt(((pi^2)/3))^(4)),null_l2=1,null_l3=0,null_l4=1/6)


xlaplace<-c(rlaplace(n=5400, location = 0, scale = 1))
NRSs(x=xlaplace,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xlaplace,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
NRSs(x=xlaplace,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=0,null_sd=sqrt(2),null_skew=0,null_kurt=(6*(sqrt(2)^4))/(4),null_l2=3/4,null_l3=0,null_l4=1/(3*sqrt(2)))


#NRSs have excellent performance even for heavy tailed distributions.

#two performance criteria, consistency (or sensitive) and standard error
library(lmom)
a=500
xpareto<-c(rpareto(5400, scale  = 1, shape=2+a/100))
targetlpareto<-lmrgpa(para = c(1,1/(2+a/100),- 1/(2+a/100)), nmom = 4)

NRSs(x=xpareto,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xpareto,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
#the standard errors are lower, especially for robust moments and L-moments
NRSs(x=xpareto,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=targetlpareto[1],null_sd=(sqrt(((2+a/100))*(1)/((-2+(2+a/100))*((-1+(2+a/100))^2)))),null_skew=((((2+a/100)+1)*(2)*(sqrt(a/100)))/((-3+(2+a/100))*(((2+a/100))^(1/2)))),null_kurt=((3+(6*((2+a/100)^3+(2+a/100)^2-6*(2+a/100)-2)/(((2+a/100))*((-3+(2+a/100)))*((-4+(2+a/100))))))),null_l2=targetlpareto[2],null_l3=targetlpareto[3],null_l4=targetlpareto[4])

a=100
xlnorm<-c(rlnorm(5400,meanlog=0,sdlog=a/100))
targetlnorm<-lmrln3(para = c(0,0, a/100), nmom = 4)

NRSs(x=xlnorm,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xlnorm,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
#the standard errors are lower, especially for robust moments and L-moments
NRSs(x=xlnorm,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=targetlnorm[1],null_sd=sqrt((exp((a/100)^2)*(-1+exp((a/100)^2)))),null_skew=sqrt(exp((a/100)^2)-1)*((2+exp((a/100)^2))),null_kurt=(((-3+exp(4*((a/100)^2))+2*exp(3*((a/100)^2))+3*exp(2*((a/100)^2))))),null_l2=targetlnorm[2],null_l3=targetlnorm[3],null_l4=targetlnorm[4])

a=150
xweibull<-c(rweibull(5400, shape=a/100, scale = 1))
library(lmom)
targetwei<-lmrwei(para = c(0, 1, a/100), nmom = 4)
NRSs(x=xweibull,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xweibull,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
NRSs(x=xweibull,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=gamma(1+1/(a/100)),null_sd=(sqrt(gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2)),null_skew=(gamma(1+3/(a/100))-3*(gamma(1+1/(a/100)))*((gamma(1+2/(a/100))))+2*((gamma(1+1/(a/100)))^3))/((sqrt(gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2))^(3)),null_kurt=((gamma(1+4/(a/100))-4*(gamma(1+3/(a/100)))*((gamma(1+1/(a/100))))+6*(gamma(1+2/(a/100)))*((gamma(1+1/(a/100)))^2)-3*((gamma(1+1/(a/100)))^4))/(((gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2))^(2))),null_l2=targetwei[2],null_l3=targetwei[3],null_l4=targetwei[4])



#for more tests, use the codes in consistency.R
