

#NRS

#I combined all the estimators into one function (easy for reviewing). There might be errors, and these are not bugs, 
#but because R is prone to producing errors for such a large function. 
#Run one function each time. If there is an error, try to restart, and then it will be fixed. 
#It will be completely fixed in the future by rewriting the code in C++.

#require foreach and doparallel for parallel processing of bootstrap (not available for some types of computers)
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
mmme<-function(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,type=1){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  etm1<-etm(sortedx,interval=interval,fast=fast,batch=batch)
  if(etm1[2]==Inf){
    return(print("ETM is infinity, due to the double precision floating point limits. Usually, the solution is transforming your original data."))
  }
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
  output1<-c(mean=mean(sortedx),etm=etm1[2],rm=rm1,qm=qm1)
  return(output1[type])
}

mmm<-function(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  etm1<-etm(sortedx,interval=interval,fast=fast,batch=batch)
  if(etm1[2]==Inf){
    return(print("ETM is infinity, due to the double precision floating point limits. Usually, the solution is transforming your original data."))
  }
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
  output1<-c(mean=mean(sortedx),etm=etm1[2],rm=rm1,qm=qm1)
  return(output1)
}
rqmean<-function (x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = FALSE,alpha = 0.05,nboot = 1000){
  if(cise){
    return (mmmci(x, interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm,alpha = alpha,nboot = nboot))
  } 
  else {return (mmm(x,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm))}
}
mmmci<-function(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,alpha=0.05,nboot=1000){
  data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  batchresults<-apply(data,1,mmm,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  bootlist1<-sort(as.matrix(batchresults[1,]))
  bootlist2<-sort(as.matrix(batchresults[2,]))
  bootlist3<-sort(as.matrix(batchresults[3,]))
  bootlist4<-sort(as.matrix(batchresults[4,]))
  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  estimate=mmm(x,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  result <- list(estimate=estimate,cimean=c(bootlist1[low],bootlist1[up]), cietm=c(bootlist2[low],bootlist2[up]),cirm=c(bootlist3[low],bootlist3[up]),
                 ciqm=c(bootlist4[low],bootlist4[up]),semean=sd(bootlist1),seetm=sd(bootlist2),serm=sd(bootlist3),seqm=sd(bootlist4))
  return(result)
}

rqscale<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.3659,drm=0.7930,dlqm=0.8218,dqm=0.7825,sd=FALSE){
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
  lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
  mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  if(sd){
    lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drm,dqm=drm)
    msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drm,dqm=drm)
    allmo<-c(l2=lmo[1],etl2=lmo[2],rl2=lmo[3],ql2=lmo[4],
           var=mo[1],etvar=mo[2],rvar=mo[3],qvar=mo[4])
    allsd<-c(l2sd=lmsd[1],etl2sd=lmsd[2],rl2sd=lmsd[3],ql2sd=lmsd[4],
             varsd=msd[1],etvarsd=msd[2],rvarsd=msd[3],qvarsd=msd[4])
    all<-c(allmo,allsd)
  }else{
    all<-c(l2=lmo[1],etl2=lmo[2],rl2=lmo[3],ql2=lmo[4],
      var=mo[1],etvar=mo[2],rvar=mo[3],qvar=mo[4])
  }
  return(all)
}

rqsd<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,drm=0.7930,dqm=0.7825){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  lengthn<-length(sortedx)
  if (boot){
    subtract<-t(replicate(times , sort(sample(sortedx, size = 2))))
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
    dp2m<-(dp^2)/2
  }
  mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  all<-c(sd=sqrt(mo[1]),etsd=sqrt(mo[2]),rsd=sqrt(mo[3]),qsd=sqrt(mo[4]))
  return(all)
}

rqtm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.1808,drm=1.7492,dlqm=1.1753,dqm=0.5715,sd=FALSE,drsd=0.7930,dqsd=0.7825){
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
  lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
  mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  if(sd){
    lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
    msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
    allmo<-c(l3=lmo[1],etl3=lmo[2],rl3=lmo[3],ql3=lmo[4],
             tm=mo[1],ettm=mo[2],rtm=mo[3],qtm=mo[4])
    allsd<-c(l3sd=lmsd[1],etl3sd=lmsd[2],rl3sd=lmsd[3],ql3sd=lmsd[4],
             tmsd=msd[1],ettmsd=msd[2],rtmsd=msd[3],qtmsd=msd[4])
    all<-c(allmo,allsd)
  }else{
    all<-c(l3=lmo[1],etl3=lmo[2],rl3=lmo[3],ql3=lmo[4],
           tm=mo[1],ettm=mo[2],rtm=mo[3],qtm=mo[4])
  }
  return(all)
}

rqfm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=-0.3542,drm=3.4560,dlqm=NaN,dqm=0.1246,sd=FALSE,drsd=0.7930,dqsd=0.7825){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  lengthn<-length(sortedx)
  
  getm<-function(vector){ 
    resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                  4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                  vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                  6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                  vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                               6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
    return(resd)
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
  
  lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
  mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  if(sd){
    lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
    msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
    allmo<-c(l4=lmo[1],etl4=lmo[2],rl4=lmo[3],ql4=lmo[4],
             fm=mo[1],etfm=mo[2],rfm=mo[3],qfm=mo[4])
    allsd<-c(l4sd=lmsd[1],etl4sd=lmsd[2],rl4sd=lmsd[3],ql4sd=lmsd[4],
             fmsd=msd[1],etfmsd=msd[2],rfmsd=msd[3],qfmsd=msd[4])
    all<-c(allmo,allsd)
  }else{
    all<-c(l4=lmo[1],etl4=lmo[2],rl4=lmo[3],ql4=lmo[4],
           fm=mo[1],etfm=mo[2],rfm=mo[3],qfm=mo[4])
  }
  return(all)
}

rLaplace<-function (n,location,scale) {
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
rPareto<-function (n, scale, shape) {
  sample1 <- scale*(runif(n))^(-1/shape)
  sample1[scale <= 0] <- NaN
  sample1[shape <= 0] <- NaN
  sample1
}

NRSssimple<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray"),sd=FALSE){
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
  mmm1<-mmm(x=sortedx,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  rqscale1<-rqscale(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=sd)
  rqtm1<-rqtm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=sd,drsd=drmscale,dqsd=dqmscale)
  rqfm1<-rqfm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=sd,drsd=drmscale,dqsd=dqmscale)
  
  if(sd){
    first<-c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4])
    second<-c(l2=rqscale1[1],etl2=rqscale1[2],rl2=rqscale1[3],ql2=rqscale1[4],sd=sqrt(rqscale1[5]),etsd=sqrt(rqscale1[6]),
              rsd=sqrt(rqscale1[7]),qsd=sqrt(rqscale1[8]))
    third<-c(l3=rqtm1[1]/rqscale1[1],etl3=rqtm1[2]/rqscale1[2],rl3=rqtm1[3]/rqscale1[3],ql3=rqtm1[4]/rqscale1[4],
             skew=(rqtm1[5])/((rqscale1[5])^(3/2)),etskew=(rqtm1[6])/((rqscale1[6])^(3/2)),
             rskew=(rqtm1[7])/((rqscale1[7])^(3/2)),qskew=(rqtm1[8])/((rqscale1[8])^(3/2)))
    fourth<-c(l4=rqfm1[1]/rqscale1[1],etl4=rqfm1[2]/rqscale1[2],rl4=rqfm1[3]/rqscale1[3],ql4=rqfm1[4]/rqscale1[4],kurt=(rqfm1[5])/((rqscale1[5])^(2)),etkurt=(rqfm1[6])/((rqscale1[6])^(2)),rkurt=(rqfm1[7])/((rqscale1[7])^(2)),qkurt=(rqfm1[8])/((rqscale1[8])^(2)))
    
    allmo<-list(first=first,second=second,third=third,fourth=fourth)
    
    
    firstsd<-c(sd=sqrt(rqscale1[5]),etsd=sqrt(rqscale1[6]),rsd=sqrt(rqscale1[7]),qsd=sqrt(rqscale1[8]))
    secondsd<-c(l2sd=rqscale1[9],etl2sd=rqscale1[10],rl2sd=rqscale1[11],ql2sd=rqscale1[12],
                sdsd=second[5]*(1/2)*(rqscale1[13]/rqscale1[5]),
                etsdsd=second[6]*(1/2)*(rqscale1[14]/rqscale1[6]),rsdsd=second[7]*(1/2)*(rqscale1[15]/rqscale1[7]),qsdsd=second[8]*(1/2)*(rqscale1[16]/rqscale1[8]))
    
    thirdsd<-c(l3sd=(rqtm1[1]/rqscale1[1])*((rqtm1[9]/rqtm1[1])^2+(rqscale1[9]/rqscale1[1])^2)^(1/2),
               etl3sd=(rqtm1[2]/rqscale1[2])*((rqtm1[10]/rqtm1[2])^2+(rqscale1[10]/rqscale1[2])^2)^(1/2),
               rl3sd=(rqtm1[3]/rqscale1[3])*((rqtm1[11]/rqtm1[3])^2+(rqscale1[11]/rqscale1[3])^2)^(1/2),
               ql3sd=(rqtm1[4]/rqscale1[4])*((rqtm1[12]/rqtm1[4])^2+(rqscale1[12]/rqscale1[4])^2)^(1/2),
               skewsd=third[5]*((rqtm1[13]/rqtm1[5])^2+(((((rqscale1[5])^(3/2))*(3/2)*(rqscale1[13]/rqscale1[5]))/((rqscale1[5])^(3/2))))^2)^(1/2),
               etskewsd=third[6]*((rqtm1[14]/rqtm1[6])^2+(((((rqscale1[6])^(3/2))*(3/2)*(rqscale1[14]/rqscale1[6]))/((rqscale1[6])^(3/2))))^2)^(1/2),
               rskewsd=third[7]*((rqtm1[15]/rqtm1[7])^2+(((((rqscale1[7])^(3/2))*(3/2)*(rqscale1[15]/rqscale1[7]))/((rqscale1[7])^(3/2))))^2)^(1/2),
               qskewsd=third[8]*((rqtm1[16]/rqtm1[8])^2+(((((rqscale1[8])^(3/2))*(3/2)*(rqscale1[16]/rqscale1[8]))/((rqscale1[8])^(3/2))))^2)^(1/2))
    
    fourthsd<-c(l4sd=(rqfm1[1]/rqscale1[1])*((rqfm1[9]/rqfm1[1])^2+(rqscale1[9]/rqscale1[1])^2)^(1/2),
               etl4sd=(rqfm1[2]/rqscale1[2])*((rqfm1[10]/rqfm1[2])^2+(rqscale1[10]/rqscale1[2])^2)^(1/2),
               rl4sd=(rqfm1[3]/rqscale1[3])*((rqfm1[11]/rqfm1[3])^2+(rqscale1[11]/rqscale1[3])^2)^(1/2),
               ql4sd=(rqfm1[4]/rqscale1[4])*((rqfm1[12]/rqfm1[4])^2+(rqscale1[12]/rqscale1[4])^2)^(1/2),
               kurtsd=fourth[5]*((rqfm1[13]/rqfm1[5])^2+(((((rqscale1[5])^(2))*(2)*(rqscale1[13]/rqscale1[5]))/((rqscale1[5])^(2))))^2)^(1/2),
               etkurtsd=fourth[6]*((rqfm1[14]/rqfm1[6])^2+(((((rqscale1[6])^(2))*(2)*(rqscale1[14]/rqscale1[6]))/((rqscale1[6])^(2))))^2)^(1/2),
               rkurtsd=fourth[7]*((rqfm1[15]/rqfm1[7])^2+(((((rqscale1[7])^(2))*(2)*(rqscale1[15]/rqscale1[7]))/((rqscale1[7])^(2))))^2)^(1/2),
               qkurtsd=fourth[8]*((rqfm1[16]/rqfm1[8])^2+(((((rqscale1[8])^(2))*(2)*(rqscale1[16]/rqscale1[8]))/((rqscale1[8])^(2))))^2)^(1/2))
               
    allsd<-list(firstsd=firstsd,secondsd=secondsd,thirdsd=thirdsd,fourthsd=fourthsd)
    
    all<-c(allmo,allsd)
  }else{
    first<-c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4])
    second<-c(l2=rqscale1[1],etl2=rqscale1[2],rl2=rqscale1[3],ql2=rqscale1[4],sd=sqrt(rqscale1[5]),etsd=sqrt(rqscale1[6]),
              rsd=sqrt(rqscale1[7]),qsd=sqrt(rqscale1[8]))
    third<-c(l3=rqtm1[1]/rqscale1[1],etl3=rqtm1[2]/rqscale1[2],rl3=rqtm1[3]/rqscale1[3],ql3=rqtm1[4]/rqscale1[4],
             skew=(rqtm1[5])/((rqscale1[5])^(3/2)),etskew=(rqtm1[6])/((rqscale1[6])^(3/2)),
             rskew=(rqtm1[7])/((rqscale1[7])^(3/2)),qskew=(rqtm1[8])/((rqscale1[8])^(3/2)))
    fourth<-c(l4=rqfm1[1]/rqscale1[1],etl4=rqfm1[2]/rqscale1[2],rl4=rqfm1[3]/rqscale1[3],ql4=rqfm1[4]/rqscale1[4],kurt=(rqfm1[5])/((rqscale1[5])^(2)),etkurt=(rqfm1[6])/((rqscale1[6])^(2)),rkurt=(rqfm1[7])/((rqscale1[7])^(2)),qkurt=(rqfm1[8])/((rqscale1[8])^(2)))
    all<-list(first=first,second=second,third=third,fourth=fourth)
  }
  
  if((rqfm1[8])/((rqscale1[8])^(2))<4.7 & standist=="exponential"){
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
  pairs23 <- data.frame()
  pairs24 <- data.frame()
  pairs25 <- data.frame()
  pairs26 <- data.frame()
  pairs27 <- data.frame()
  pairs28 <- data.frame()
  for (i in 1:nboot) {
    mmm1<-mmm(x=data[i,],interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
    rqscale1<-rqscale(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=FALSE)
    rqtm1<-rqtm(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
    rqfm1<-rqfm(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
    estimate1<-c(c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4]),
                 c(l2=rqscale1[1],etl2=rqscale1[2],rl2=rqscale1[3],ql2=rqscale1[4],sd=sqrt(rqscale1[5]),etsd=sqrt(rqscale1[6]),
                   rsd=sqrt(rqscale1[7]),qsd=sqrt(rqscale1[8])),
                 c(l3=rqtm1[1]/rqscale1[1],etl3=rqtm1[2]/rqscale1[2],rl3=rqtm1[3]/rqscale1[3],ql3=rqtm1[4]/rqscale1[4],
                   skew=(rqtm1[5])/((rqscale1[5])^(3/2)),etskew=(rqtm1[6])/((rqscale1[6])^(3/2)),
                   rskew=(rqtm1[7])/((rqscale1[7])^(3/2)),qskew=(rqtm1[8])/((rqscale1[8])^(3/2))),
                 c(l4=rqfm1[1]/rqscale1[1],etl4=rqfm1[2]/rqscale1[2],rl4=rqfm1[3]/rqscale1[3],ql4=rqfm1[4]/rqscale1[4],kurt=(rqfm1[5])/((rqscale1[5])^(2)),etkurt=(rqfm1[6])/((rqscale1[6])^(2)),rkurt=(rqfm1[7])/((rqscale1[7])^(2)),qkurt=(rqfm1[8])/((rqscale1[8])^(2))))
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
    pairs23 <- rbind(pairs23,estimate1[23])
    pairs24 <- rbind(pairs24,estimate1[24])
    pairs25 <- rbind(pairs25,estimate1[25])
    pairs26 <- rbind(pairs26,estimate1[26])
    pairs27 <- rbind(pairs27,estimate1[27])
    pairs28 <- rbind(pairs28,estimate1[28])
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
  bootlist23<-sort(as.matrix(pairs23))
  bootlist24<-sort(as.matrix(pairs24))
  bootlist25<-sort(as.matrix(pairs25))
  bootlist26<-sort(as.matrix(pairs26))
  bootlist27<-sort(as.matrix(pairs27))
  bootlist28<-sort(as.matrix(pairs28))
  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  pbm<-function(bootlist,null_value){
    p <- mean(bootlist > null_value) + 0.5 * mean(bootlist == null_value)
    p <- 2 * min(c(p, 1 - p))
    return(p)
  }
  mmm1<-mmm(x=sortedx,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  rqscale1<-rqscale(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=FALSE)
  rqtm1<-rqtm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
  rqfm1<-rqfm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
  estimate<-c(c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4]),
               c(l2=rqscale1[1],etl2=rqscale1[2],rl2=rqscale1[3],ql2=rqscale1[4],sd=sqrt(rqscale1[5]),etsd=sqrt(rqscale1[6]),
                 rsd=sqrt(rqscale1[7]),qsd=sqrt(rqscale1[8])),
               c(l3=rqtm1[1]/rqscale1[1],etl3=rqtm1[2]/rqscale1[2],rl3=rqtm1[3]/rqscale1[3],ql3=rqtm1[4]/rqscale1[4],
                 skew=(rqtm1[5])/((rqscale1[5])^(3/2)),etskew=(rqtm1[6])/((rqscale1[6])^(3/2)),
                 rskew=(rqtm1[7])/((rqscale1[7])^(3/2)),qskew=(rqtm1[8])/((rqscale1[8])^(3/2))),
               c(l4=rqfm1[1]/rqscale1[1],etl4=rqfm1[2]/rqscale1[2],rl4=rqfm1[3]/rqscale1[3],ql4=rqfm1[4]/rqscale1[4],kurt=(rqfm1[5])/((rqscale1[5])^(2)),etkurt=(rqfm1[6])/((rqscale1[6])^(2)),rkurt=(rqfm1[7])/((rqscale1[7])^(2)),qkurt=(rqfm1[8])/((rqscale1[8])^(2))))
  p_value<-c(mean=pbm(bootlist1,null_mean),etm=pbm(bootlist2,null_mean),rm=pbm(bootlist3,null_mean),qm=pbm(bootlist4,null_mean),
             l2=pbm(bootlist5,null_l2),etl2=pbm(bootlist6,null_l2),
             rl2=pbm(bootlist7,null_l2),ql2=pbm(bootlist8,null_l2),
             sd=pbm(bootlist9,null_sd),etsd=pbm(bootlist10,null_sd),
             rsd=pbm(bootlist11,null_sd),qsd=pbm(bootlist12,null_sd),
             l3=pbm(bootlist13,null_l3),etl3=pbm(bootlist14,null_l3),rl3=pbm(bootlist15,null_l3),ql3=pbm(bootlist16,null_l3),
             
             skew=pbm(bootlist17,null_skew),etskew=pbm(bootlist18,null_skew),rskew=pbm(bootlist19,null_skew),qskew=pbm(bootlist20,null_skew),
             
             l4=pbm(bootlist21,null_l4),etl4=pbm(bootlist22,null_l4),rl4=pbm(bootlist23,null_l4),ql4=pbm(bootlist24,null_l4),
             
             kurt=pbm(bootlist25,null_kurt),etkurt=pbm(bootlist26,null_kurt),rkurt=pbm(bootlist27,null_kurt),qkurt=pbm(bootlist28,null_kurt))
  
  ci<-c(mean=c(bootlist1[low],bootlist1[up]),etm=c(bootlist2[low],bootlist2[up]),rm=c(bootlist3[low],bootlist3[up]),qm=c(bootlist4[low],bootlist4[up]),
        l2=c(bootlist5[low],bootlist5[up]),etl2=c(bootlist6[low],bootlist6[up]),rl2=c(bootlist7[low],bootlist7[up]),ql2=c(bootlist8[low],bootlist8[up]),
        
        sd=c(bootlist9[low],bootlist9[up]),etsd=c(bootlist10[low],bootlist10[up]),rsd=c(bootlist11[low],bootlist11[up]),qsd=c(bootlist12[low],bootlist12[up]),
        
        l3=c(bootlist13[low],bootlist13[up]),etl3=c(bootlist14[low],bootlist14[up]),rl3=c(bootlist15[low],bootlist15[up]),ql3=c(bootlist16[low],bootlist16[up]),
        
        skew=c(bootlist17[low],bootlist17[up]),etskew=c(bootlist18[low],bootlist18[up]),rskew=c(bootlist19[low],bootlist19[up]),qskew=c(bootlist20[low],bootlist20[up]),
        
        l4=c(bootlist21[low],bootlist21[up]),etl4=c(bootlist22[low],bootlist22[up]),rl4=c(bootlist23[low],bootlist23[up]),ql4=c(bootlist24[low],bootlist24[up]),
        
        kurt=c(bootlist25[low],bootlist25[up]),etkurt=c(bootlist26[low],bootlist26[up]),rkurt=c(bootlist27[low],bootlist27[up]),qkurt=c(bootlist28[low],bootlist28[up]))
  
  se<-c(mean=sd(bootlist1),etm=sd(bootlist2),rm=sd(bootlist3),qm=sd(bootlist4),
    l2=sd(bootlist5),etl2=sd(bootlist6),
    rl2=sd(bootlist7),ql2=sd(bootlist8),
    sd=sd(bootlist9),etsd=sd(bootlist10),
    rsd=sd(bootlist11),qsd=sd(bootlist12),
    l3=sd(bootlist13),etl3=sd(bootlist14),rl3=sd(bootlist15),ql3=sd(bootlist16),
    
    skew=sd(bootlist17),etskew=sd(bootlist18,null_skew),rskew=sd(bootlist19),qskew=sd(bootlist20),
    
    l4=sd(bootlist21),etl4=sd(bootlist22),rl4=sd(bootlist23),ql4=sd(bootlist24),
    
    kurt=sd(bootlist25),etkurt=sd(bootlist26),rkurt=sd(bootlist27),qkurt=sd(bootlist28))

  all<-list(estimate=estimate,ci=ci,se=se,p_value)
  if((rqfm1[8])/((rqscale1[8])^(2))<4.7 & standist=="exponential"){
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
    mmm<-function(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      etm1<-etm(sortedx,interval=interval,fast=fast,batch=batch)
      if(etm1[2]==Inf){
        return(print("ETM is infinity, due to the double precision floating point limits. Usually, the solution is transforming your original data."))
      }
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
      output1<-c(mean=mean(sortedx),etm=etm1[2],rm=rm1,qm=qm1)
      return(output1)
    }
    rqmean<-function (x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = FALSE,alpha = 0.05,nboot = 1000){
      if(cise){
        return (mmmci(x, interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm,alpha = alpha,nboot = nboot))
      } 
      else {return (mmm(x,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm))}
    }

    rqscale<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.3659,drm=0.7930,dlqm=0.8218,dqm=0.7825,sd=FALSE){
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
      lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      if(sd){
        lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drm,dqm=drm)
        msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drm,dqm=drm)
        allmo<-c(l2=lmo[1],etl2=lmo[2],rl2=lmo[3],ql2=lmo[4],
                 var=mo[1],etvar=mo[2],rvar=mo[3],qvar=mo[4])
        allsd<-c(l2sd=lmsd[1],etl2sd=lmsd[2],rl2sd=lmsd[3],ql2sd=lmsd[4],
                 varsd=msd[1],etvarsd=msd[2],rvarsd=msd[3],qvarsd=msd[4])
        all<-c(allmo,allsd)
      }else{
        all<-c(l2=lmo[1],etl2=lmo[2],rl2=lmo[3],ql2=lmo[4],
               var=mo[1],etvar=mo[2],rvar=mo[3],qvar=mo[4])
      }
      return(all)
    }
    
    rqsd<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,drm=0.7930,dqm=0.7825){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthn<-length(sortedx)
      if (boot){
        subtract<-t(replicate(times , sort(sample(sortedx, size = 2))))
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
        dp2m<-(dp^2)/2
      }
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      all<-c(sd=sqrt(mo[1]),etsd=sqrt(mo[2]),rsd=sqrt(mo[3]),qsd=sqrt(mo[4]))
      return(all)
    }
    
    rqtm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.1808,drm=1.7492,dlqm=1.1753,dqm=0.5715,sd=FALSE,drsd=0.7930,dqsd=0.7825){
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
      lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      if(sd){
        lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        allmo<-c(l3=lmo[1],etl3=lmo[2],rl3=lmo[3],ql3=lmo[4],
                 tm=mo[1],ettm=mo[2],rtm=mo[3],qtm=mo[4])
        allsd<-c(l3sd=lmsd[1],etl3sd=lmsd[2],rl3sd=lmsd[3],ql3sd=lmsd[4],
                 tmsd=msd[1],ettmsd=msd[2],rtmsd=msd[3],qtmsd=msd[4])
        all<-c(allmo,allsd)
      }else{
        all<-c(l3=lmo[1],etl3=lmo[2],rl3=lmo[3],ql3=lmo[4],
               tm=mo[1],ettm=mo[2],rtm=mo[3],qtm=mo[4])
      }
      return(all)
    }
    
    rqfm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=-0.3542,drm=3.4560,dlqm=NaN,dqm=0.1246,sd=FALSE,drsd=0.7930,dqsd=0.7825){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthn<-length(sortedx)
      
      getm<-function(vector){ 
        resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                      4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                      vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                      6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                      vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                                   6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
        return(resd)
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
      
      lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      if(sd){
        lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        allmo<-c(l4=lmo[1],etl4=lmo[2],rl4=lmo[3],ql4=lmo[4],
                 fm=mo[1],etfm=mo[2],rfm=mo[3],qfm=mo[4])
        allsd<-c(l4sd=lmsd[1],etl4sd=lmsd[2],rl4sd=lmsd[3],ql4sd=lmsd[4],
                 fmsd=msd[1],etfmsd=msd[2],rfmsd=msd[3],qfmsd=msd[4])
        all<-c(allmo,allsd)
      }else{
        all<-c(l4=lmo[1],etl4=lmo[2],rl4=lmo[3],ql4=lmo[4],
               fm=mo[1],etfm=mo[2],rfm=mo[3],qfm=mo[4])
      }
      return(all)
    }
    mmm1<-mmm(x=data[i,],interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
    rqscale1<-rqscale(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=FALSE)
    rqtm1<-rqtm(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
    rqfm1<-rqfm(x=data[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
    estimate1<-c(c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4]),
                 c(l2=rqscale1[1],etl2=rqscale1[2],rl2=rqscale1[3],ql2=rqscale1[4],sd=sqrt(rqscale1[5]),etsd=sqrt(rqscale1[6]),
                   rsd=sqrt(rqscale1[7]),qsd=sqrt(rqscale1[8])),
                 c(l3=rqtm1[1]/rqscale1[1],etl3=rqtm1[2]/rqscale1[2],rl3=rqtm1[3]/rqscale1[3],ql3=rqtm1[4]/rqscale1[4],
                   skew=(rqtm1[5])/((rqscale1[5])^(3/2)),etskew=(rqtm1[6])/((rqscale1[6])^(3/2)),
                   rskew=(rqtm1[7])/((rqscale1[7])^(3/2)),qskew=(rqtm1[8])/((rqscale1[8])^(3/2))),
                 c(l4=rqfm1[1]/rqscale1[1],etl4=rqfm1[2]/rqscale1[2],rl4=rqfm1[3]/rqscale1[3],ql4=rqfm1[4]/rqscale1[4],kurt=(rqfm1[5])/((rqscale1[5])^(2)),etkurt=(rqfm1[6])/((rqscale1[6])^(2)),rkurt=(rqfm1[7])/((rqscale1[7])^(2)),qkurt=(rqfm1[8])/((rqscale1[8])^(2))))
    
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
  bootlist23<-sort(as.matrix(estimate0[,23]))
  bootlist24<-sort(as.matrix(estimate0[,24]))
  bootlist25<-sort(as.matrix(estimate0[,25]))
  bootlist26<-sort(as.matrix(estimate0[,26]))
  bootlist27<-sort(as.matrix(estimate0[,27]))
  bootlist28<-sort(as.matrix(estimate0[,28]))
  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  pbm<-function(bootlist,null_value){
    p <- mean(bootlist > null_value) + 0.5 * mean(bootlist == null_value)
    p <- 2 * min(c(p, 1 - p))
    return(p)
  }
  mmm1<-mmm(x=sortedx,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  rqscale1<-rqscale(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=FALSE)
  rqtm1<-rqtm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
  rqfm1<-rqfm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
  estimate<-c(c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4]),
              c(l2=rqscale1[1],etl2=rqscale1[2],rl2=rqscale1[3],ql2=rqscale1[4],sd=sqrt(rqscale1[5]),etsd=sqrt(rqscale1[6]),
                rsd=sqrt(rqscale1[7]),qsd=sqrt(rqscale1[8])),
              c(l3=rqtm1[1]/rqscale1[1],etl3=rqtm1[2]/rqscale1[2],rl3=rqtm1[3]/rqscale1[3],ql3=rqtm1[4]/rqscale1[4],
                skew=(rqtm1[5])/((rqscale1[5])^(3/2)),etskew=(rqtm1[6])/((rqscale1[6])^(3/2)),
                rskew=(rqtm1[7])/((rqscale1[7])^(3/2)),qskew=(rqtm1[8])/((rqscale1[8])^(3/2))),
              c(l4=rqfm1[1]/rqscale1[1],etl4=rqfm1[2]/rqscale1[2],rl4=rqfm1[3]/rqscale1[3],ql4=rqfm1[4]/rqscale1[4],kurt=(rqfm1[5])/((rqscale1[5])^(2)),etkurt=(rqfm1[6])/((rqscale1[6])^(2)),rkurt=(rqfm1[7])/((rqscale1[7])^(2)),qkurt=(rqfm1[8])/((rqscale1[8])^(2))))
  p_value<-c(mean=pbm(bootlist1,null_mean),etm=pbm(bootlist2,null_mean),rm=pbm(bootlist3,null_mean),qm=pbm(bootlist4,null_mean),
             l2=pbm(bootlist5,null_l2),etl2=pbm(bootlist6,null_l2),
             rl2=pbm(bootlist7,null_l2),ql2=pbm(bootlist8,null_l2),
             sd=pbm(bootlist9,null_sd),etsd=pbm(bootlist10,null_sd),
             rsd=pbm(bootlist11,null_sd),qsd=pbm(bootlist12,null_sd),
             l3=pbm(bootlist13,null_l3),etl3=pbm(bootlist14,null_l3),rl3=pbm(bootlist15,null_l3),ql3=pbm(bootlist16,null_l3),
             
             skew=pbm(bootlist17,null_skew),etskew=pbm(bootlist18,null_skew),rskew=pbm(bootlist19,null_skew),qskew=pbm(bootlist20,null_skew),
             
             l4=pbm(bootlist21,null_l4),etl4=pbm(bootlist22,null_l4),rl4=pbm(bootlist23,null_l4),ql4=pbm(bootlist24,null_l4),
             
             kurt=pbm(bootlist25,null_kurt),etkurt=pbm(bootlist26,null_kurt),rkurt=pbm(bootlist27,null_kurt),qkurt=pbm(bootlist28,null_kurt))
  
  ci<-c(mean=c(bootlist1[low],bootlist1[up]),etm=c(bootlist2[low],bootlist2[up]),rm=c(bootlist3[low],bootlist3[up]),qm=c(bootlist4[low],bootlist4[up]),
        l2=c(bootlist5[low],bootlist5[up]),etl2=c(bootlist6[low],bootlist6[up]),rl2=c(bootlist7[low],bootlist7[up]),ql2=c(bootlist8[low],bootlist8[up]),
        
        sd=c(bootlist9[low],bootlist9[up]),etsd=c(bootlist10[low],bootlist10[up]),rsd=c(bootlist11[low],bootlist11[up]),qsd=c(bootlist12[low],bootlist12[up]),
        
        l3=c(bootlist13[low],bootlist13[up]),etl3=c(bootlist14[low],bootlist14[up]),rl3=c(bootlist15[low],bootlist15[up]),ql3=c(bootlist16[low],bootlist16[up]),
        
        skew=c(bootlist17[low],bootlist17[up]),etskew=c(bootlist18[low],bootlist18[up]),rskew=c(bootlist19[low],bootlist19[up]),qskew=c(bootlist20[low],bootlist20[up]),
        
        l4=c(bootlist21[low],bootlist21[up]),etl4=c(bootlist22[low],bootlist22[up]),rl4=c(bootlist23[low],bootlist23[up]),ql4=c(bootlist24[low],bootlist24[up]),
        
        kurt=c(bootlist25[low],bootlist25[up]),etkurt=c(bootlist26[low],bootlist26[up]),rkurt=c(bootlist27[low],bootlist27[up]),qkurt=c(bootlist28[low],bootlist28[up]))
  
  se<-c(mean=sd(bootlist1),etm=sd(bootlist2),rm=sd(bootlist3),qm=sd(bootlist4),
        l2=sd(bootlist5),etl2=sd(bootlist6),
        rl2=sd(bootlist7),ql2=sd(bootlist8),
        sd=sd(bootlist9),etsd=sd(bootlist10),
        rsd=sd(bootlist11),qsd=sd(bootlist12),
        l3=sd(bootlist13),etl3=sd(bootlist14),rl3=sd(bootlist15),ql3=sd(bootlist16),
        
        skew=sd(bootlist17),etskew=sd(bootlist18,null_skew),rskew=sd(bootlist19),qskew=sd(bootlist20),
        
        l4=sd(bootlist21),etl4=sd(bootlist22),rl4=sd(bootlist23),ql4=sd(bootlist24),
        
        kurt=sd(bootlist25),etkurt=sd(bootlist26),rkurt=sd(bootlist27),qkurt=sd(bootlist28))
  
  all<-list(estimate=estimate,ci=ci,se=se,p_value)
  if((rqfm1[8])/((rqscale1[8])^(2))<4.7 & standist=="exponential"){
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
    mmm<-function(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      etm1<-etm(sortedx,interval=interval,fast=fast,batch=batch)
      if(etm1[2]==Inf){
        return(print("ETM is infinity, due to the double precision floating point limits. Usually, the solution is transforming your original data."))
      }
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
      output1<-c(mean=mean(sortedx),etm=etm1[2],rm=rm1,qm=qm1)
      return(output1)
    }
    rqmean<-function (x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = FALSE,alpha = 0.05,nboot = 1000){
      if(cise){
        return (mmmci(x, interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm,alpha = alpha,nboot = nboot))
      } 
      else {return (mmm(x,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm))}
    }
    
    rqscale<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.3659,drm=0.7930,dlqm=0.8218,dqm=0.7825,sd=FALSE){
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
      lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      if(sd){
        lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drm,dqm=drm)
        msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drm,dqm=drm)
        allmo<-c(l2=lmo[1],etl2=lmo[2],rl2=lmo[3],ql2=lmo[4],
                 var=mo[1],etvar=mo[2],rvar=mo[3],qvar=mo[4])
        allsd<-c(l2sd=lmsd[1],etl2sd=lmsd[2],rl2sd=lmsd[3],ql2sd=lmsd[4],
                 varsd=msd[1],etvarsd=msd[2],rvarsd=msd[3],qvarsd=msd[4])
        all<-c(allmo,allsd)
      }else{
        all<-c(l2=lmo[1],etl2=lmo[2],rl2=lmo[3],ql2=lmo[4],
               var=mo[1],etvar=mo[2],rvar=mo[3],qvar=mo[4])
      }
      return(all)
    }
    
    rqsd<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,drm=0.7930,dqm=0.7825){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthn<-length(sortedx)
      if (boot){
        subtract<-t(replicate(times , sort(sample(sortedx, size = 2))))
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
        dp2m<-(dp^2)/2
      }
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      all<-c(sd=sqrt(mo[1]),etsd=sqrt(mo[2]),rsd=sqrt(mo[3]),qsd=sqrt(mo[4]))
      return(all)
    }
    
    rqtm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.1808,drm=1.7492,dlqm=1.1753,dqm=0.5715,sd=FALSE,drsd=0.7930,dqsd=0.7825){
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
      lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      if(sd){
        lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        allmo<-c(l3=lmo[1],etl3=lmo[2],rl3=lmo[3],ql3=lmo[4],
                 tm=mo[1],ettm=mo[2],rtm=mo[3],qtm=mo[4])
        allsd<-c(l3sd=lmsd[1],etl3sd=lmsd[2],rl3sd=lmsd[3],ql3sd=lmsd[4],
                 tmsd=msd[1],ettmsd=msd[2],rtmsd=msd[3],qtmsd=msd[4])
        all<-c(allmo,allsd)
      }else{
        all<-c(l3=lmo[1],etl3=lmo[2],rl3=lmo[3],ql3=lmo[4],
               tm=mo[1],ettm=mo[2],rtm=mo[3],qtm=mo[4])
      }
      return(all)
    }
    
    rqfm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=-0.3542,drm=3.4560,dlqm=NaN,dqm=0.1246,sd=FALSE,drsd=0.7930,dqsd=0.7825){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthn<-length(sortedx)
      
      getm<-function(vector){ 
        resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                      4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                      vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                      6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                      vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                                   6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
        return(resd)
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
      
      lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      if(sd){
        lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        allmo<-c(l4=lmo[1],etl4=lmo[2],rl4=lmo[3],ql4=lmo[4],
                 fm=mo[1],etfm=mo[2],rfm=mo[3],qfm=mo[4])
        allsd<-c(l4sd=lmsd[1],etl4sd=lmsd[2],rl4sd=lmsd[3],ql4sd=lmsd[4],
                 fmsd=msd[1],etfmsd=msd[2],rfmsd=msd[3],qfmsd=msd[4])
        all<-c(allmo,allsd)
      }else{
        all<-c(l4=lmo[1],etl4=lmo[2],rl4=lmo[3],ql4=lmo[4],
               fm=mo[1],etfm=mo[2],rfm=mo[3],qfm=mo[4])
      }
      return(all)
    }
    mmmx<-mmm(x=datax[i,],interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
    rqscalex<-rqscale(x=datax[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=FALSE)
    rqtmx<-rqtm(x=datax[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
    rqfmx<-rqfm(x=datax[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
    mmmy<-mmm(x=datay[i,],interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
    rqscaley<-rqscale(x=datay[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=FALSE)
    rqtmy<-rqtm(x=datay[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
    rqfmy<-rqfm(x=datay[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
    
    estimate1<-c(c(meanx=mmmx[1],etmx=mmmx[2],rmx=mmmx[3],qmx=mmmx[4]),
                 c(l2x=rqscalex[1],etl2x=rqscalex[2],rl2x=rqscalex[3],ql2x=rqscalex[4],sdx=sqrt(rqscalex[5]),etsdx=sqrt(rqscalex[6]),
                   rsdx=sqrt(rqscalex[7]),qsdx=sqrt(rqscalex[8])),
                 c(l3x=rqtmx[1]/rqscalex[1],etl3x=rqtmx[2]/rqscalex[2],rl3x=rqtmx[3]/rqscalex[3],ql3x=rqtmx[4]/rqscalex[4],
                   skewx=(rqtmx[5])/((rqscalex[5])^(3/2)),etskewx=(rqtmx[6])/((rqscalex[6])^(3/2)),
                   rskewx=(rqtmx[7])/((rqscalex[7])^(3/2)),qskewx=(rqtmx[8])/((rqscalex[8])^(3/2))),
                 c(l4x=rqfmx[1]/rqscalex[1],etl4x=rqfmx[2]/rqscalex[2],rl4x=rqfmx[3]/rqscalex[3],ql4x=rqfmx[4]/rqscalex[4],kurtx=(rqfmx[5])/((rqscalex[5])^(2)),etkurtx=(rqfmx[6])/((rqscalex[6])^(2)),rkurtx=(rqfmx[7])/((rqscalex[7])^(2)),qkurtx=(rqfmx[8])/((rqscalex[8])^(2))),
                 c(meany=mmmy[1],etmy=mmmy[2],rmy=mmmy[3],qmy=mmmy[4]),
                 c(l2y=rqscaley[1],etl2y=rqscaley[2],rl2y=rqscaley[3],ql2y=rqscaley[4],sdy=sqrt(rqscaley[5]),etsdy=sqrt(rqscaley[6]),
                   rsdy=sqrt(rqscaley[7]),qsdy=sqrt(rqscaley[8])),
                 c(l3y=rqtmy[1]/rqscaley[1],etl3y=rqtmy[2]/rqscaley[2],rl3y=rqtmy[3]/rqscaley[3],ql3y=rqtmy[4]/rqscaley[4],
                   skewy=(rqtmy[5])/((rqscaley[5])^(3/2)),etskewy=(rqtmy[6])/((rqscaley[6])^(3/2)),
                   rskewy=(rqtmy[7])/((rqscaley[7])^(3/2)),qskewy=(rqtmy[8])/((rqscaley[8])^(3/2))),
                 c(l4y=rqfmy[1]/rqscaley[1],etl4y=rqfmy[2]/rqscaley[2],rl4y=rqfmy[3]/rqscaley[3],ql4y=rqfmy[4]/rqscaley[4],kurty=(rqfmy[5])/((rqscaley[5])^(2)),etkurty=(rqfmy[6])/((rqscaley[6])^(2)),rkurty=(rqfmy[7])/((rqscaley[7])^(2)),qkurty=(rqfmy[8])/((rqscaley[8])^(2)))
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
  bootlist45<-(as.matrix(estimate0[,45]))
  bootlist46<-(as.matrix(estimate0[,46]))
  bootlist47<-(as.matrix(estimate0[,47]))
  bootlist48<-(as.matrix(estimate0[,48]))
  bootlist49<-(as.matrix(estimate0[,49]))
  bootlist50<-(as.matrix(estimate0[,50]))
  bootlist51<-(as.matrix(estimate0[,51]))
  bootlist52<-(as.matrix(estimate0[,52]))
  bootlist53<-(as.matrix(estimate0[,53]))
  bootlist54<-(as.matrix(estimate0[,54]))
  bootlist55<-(as.matrix(estimate0[,55]))
  bootlist56<-(as.matrix(estimate0[,56]))

  bootlist1a<-sort(bootlist1-bootlist29)
  bootlist2a<-sort(bootlist2-bootlist30)
  bootlist3a<-sort(bootlist3-bootlist31)
  bootlist4a<-sort(bootlist4-bootlist32)
  bootlist5a<-sort(bootlist5-bootlist33)
  bootlist6a<-sort(bootlist6-bootlist34)
  bootlist7a<-sort(bootlist7-bootlist35)
  bootlist8a<-sort(bootlist8-bootlist36)
  bootlist9a<-sort(bootlist9-bootlist37)
  bootlist10a<-sort(bootlist10-bootlist38)
  bootlist11a<-sort(bootlist11-bootlist39)
  bootlist12a<-sort(bootlist12-bootlist40)
  bootlist13a<-sort(bootlist13-bootlist41)
  bootlist14a<-sort(bootlist14-bootlist42)
  bootlist15a<-sort(bootlist15-bootlist43)
  bootlist16a<-sort(bootlist16-bootlist44)
  bootlist17a<-sort(bootlist17-bootlist45)
  bootlist18a<-sort(bootlist18-bootlist46)
  bootlist19a<-sort(bootlist19-bootlist47)
  bootlist20a<-sort(bootlist20-bootlist48)
  bootlist21a<-sort(bootlist21-bootlist49)
  bootlist22a<-sort(bootlist22-bootlist50)
  bootlist23a<-sort(bootlist23-bootlist51)
  bootlist24a<-sort(bootlist24-bootlist52)
  bootlist25a<-sort(bootlist25-bootlist53)
  bootlist26a<-sort(bootlist26-bootlist54)
  bootlist27a<-sort(bootlist27-bootlist55)
  bootlist28a<-sort(bootlist28-bootlist56)
  
  
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
  bootlist45<-sort(bootlist45)
  bootlist46<-sort(bootlist46)
  bootlist47<-sort(bootlist47)
  bootlist48<-sort(bootlist48)
  bootlist49<-sort(bootlist49)
  bootlist50<-sort(bootlist50)
  bootlist51<-sort(bootlist51)
  bootlist52<-sort(bootlist52)
  bootlist53<-sort(bootlist53)
  bootlist54<-sort(bootlist54)
  bootlist55<-sort(bootlist55)
  bootlist56<-sort(bootlist56)

  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  
  pb2<-function(bootlist){
    p <- mean(bootlist < 0) + 0.5 * mean(bootlist == 0)
    p <- 2 * min(c(p, 1 - p))
    p
  }
  cidiff<-function(bootlist,low,up){
    c(bootlist[low],bootlist[up])
  }
  ci_diff<-c(c(mean=cidiff(bootlist1a,low=low,up=up),etm=cidiff(bootlist2a,low=low,up=up),rm=cidiff(bootlist3a,low=low,up=up),qm=cidiff(bootlist4a,low=low,up=up)),
             c(l2=cidiff(bootlist5a,low=low,up=up),etl2=cidiff(bootlist6a,low=low,up=up),rl2=cidiff(bootlist7a,low=low,up=up),ql2=cidiff(bootlist8a,low=low,up=up),sd=cidiff(bootlist9a,low=low,up=up),etsd=cidiff(bootlist10a,low=low,up=up),
               rsd=cidiff(bootlist11a,low=low,up=up),qsd=cidiff(bootlist12a,low=low,up=up)),
             c(l3=cidiff(bootlist13a,low=low,up=up),etl3=cidiff(bootlist14a,low=low,up=up),rl3=cidiff(bootlist15a,low=low,up=up),ql3=cidiff(bootlist16a,low=low,up=up),
               skew=cidiff(bootlist17a,low=low,up=up),etskew=cidiff(bootlist18a,low=low,up=up),
               rskew=cidiff(bootlist19a,low=low,up=up),qskew=cidiff(bootlist20a,low=low,up=up)),
             c(l4=cidiff(bootlist21a,low=low,up=up),etl4=cidiff(bootlist22a,low=low,up=up),rl4=cidiff(bootlist23a,low=low,up=up),ql4=cidiff(bootlist24a,low=low,up=up),kurt=cidiff(bootlist25a,low=low,up=up),etkurt=cidiff(bootlist26a,low=low,up=up),rkurt=cidiff(bootlist27a,low=low,up=up),qkurt=cidiff(bootlist28a,low=low,up=up)))
    
  
  p_value_diff<-c(c(mean=pb2(bootlist1a),etm=pb2(bootlist2a),rm=pb2(bootlist3a),qm=pb2(bootlist4a)),
                c(l2=pb2(bootlist5a),etl2=pb2(bootlist6a),rl2=pb2(bootlist7a),ql2=pb2(bootlist8a),sd=pb2(bootlist9a),etsd=pb2(bootlist10a),
                  rsd=pb2(bootlist11a),qsd=pb2(bootlist12a)),
                c(l3=pb2(bootlist13a),etl3=pb2(bootlist14a),rl3=pb2(bootlist15a),ql3=pb2(bootlist16a),
                  skew=pb2(bootlist17a),etskew=pb2(bootlist18a),
                  rskew=pb2(bootlist19a),qskew=pb2(bootlist20a)),
                c(l4=pb2(bootlist21a),etl4=pb2(bootlist22a),rl4=pb2(bootlist23a),ql4=pb2(bootlist24a),kurt=pb2(bootlist25a),etkurt=pb2(bootlist26a),rkurt=pb2(bootlist27a),qkurt=pb2(bootlist28a)))
  
  
  mmmx<-mmm(x=sortedx,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  rqscalex<-rqscale(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=FALSE)
  rqtmx<-rqtm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
  rqfmx<-rqfm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
  mmmy<-mmm(x=sortedy,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  rqscaley<-rqscale(x=sortedy,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=FALSE)
  rqtmy<-rqtm(x=sortedy,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
  rqfmy<-rqfm(x=sortedy,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
  
  estimate1<-c(c(meanx=mmmx[1],etmx=mmmx[2],rmx=mmmx[3],qmx=mmmx[4]),
               c(l2x=rqscalex[1],etl2x=rqscalex[2],rl2x=rqscalex[3],ql2x=rqscalex[4],sdx=sqrt(rqscalex[5]),etsdx=sqrt(rqscalex[6]),
                 rsdx=sqrt(rqscalex[7]),qsdx=sqrt(rqscalex[8])),
               c(l3x=rqtmx[1]/rqscalex[1],etl3x=rqtmx[2]/rqscalex[2],rl3x=rqtmx[3]/rqscalex[3],ql3x=rqtmx[4]/rqscalex[4],
                 skewx=(rqtmx[5])/((rqscalex[5])^(3/2)),etskewx=(rqtmx[6])/((rqscalex[6])^(3/2)),
                 rskewx=(rqtmx[7])/((rqscalex[7])^(3/2)),qskewx=(rqtmx[8])/((rqscalex[8])^(3/2))),
               c(l4x=rqfmx[1]/rqscalex[1],etl4x=rqfmx[2]/rqscalex[2],rl4x=rqfmx[3]/rqscalex[3],ql4x=rqfmx[4]/rqscalex[4],kurtx=(rqfmx[5])/((rqscalex[5])^(2)),etkurtx=(rqfmx[6])/((rqscalex[6])^(2)),rkurtx=(rqfmx[7])/((rqscalex[7])^(2)),qkurtx=(rqfmx[8])/((rqscalex[8])^(2))),
               c(meany=mmmy[1],etmy=mmmy[2],rmy=mmmy[3],qmy=mmmy[4]),
               c(l2y=rqscaley[1],etl2y=rqscaley[2],rl2y=rqscaley[3],ql2y=rqscaley[4],sdy=sqrt(rqscaley[5]),etsdy=sqrt(rqscaley[6]),
                 rsdy=sqrt(rqscaley[7]),qsdy=sqrt(rqscaley[8])),
               c(l3y=rqtmy[1]/rqscaley[1],etl3y=rqtmy[2]/rqscaley[2],rl3y=rqtmy[3]/rqscaley[3],ql3y=rqtmy[4]/rqscaley[4],
                 skewy=(rqtmy[5])/((rqscaley[5])^(3/2)),etskewy=(rqtmy[6])/((rqscaley[6])^(3/2)),
                 rskewy=(rqtmy[7])/((rqscaley[7])^(3/2)),qskewy=(rqtmy[8])/((rqscaley[8])^(3/2))),
               c(l4y=rqfmy[1]/rqscaley[1],etl4y=rqfmy[2]/rqscaley[2],rl4y=rqfmy[3]/rqscaley[3],ql4y=rqfmy[4]/rqscaley[4],kurty=(rqfmy[5])/((rqscaley[5])^(2)),etkurty=(rqfmy[6])/((rqscaley[6])^(2)),rkurty=(rqfmy[7])/((rqscaley[7])^(2)),qkurty=(rqfmy[8])/((rqscaley[8])^(2)))
  )
  
  ci<-c(meanx=c(bootlist1[low],bootlist1[up]),etmx=c(bootlist2[low],bootlist2[up]),rmx=c(bootlist3[low],bootlist3[up]),qmx=c(bootlist4[low],bootlist4[up]),
        l2x=c(bootlist5[low],bootlist5[up]),etl2x=c(bootlist6[low],bootlist6[up]),rl2x=c(bootlist7[low],bootlist7[up]),ql2x=c(bootlist8[low],bootlist8[up]),
        sdx=c(bootlist9[low],bootlist9[up]),
        etsdx=c(bootlist10[low],bootlist10[up]),
        rsdx=c(bootlist11[low],bootlist11[up]),qsdx=c(bootlist12[low],bootlist12[up]),
        l3x=c(bootlist13[low],bootlist13[up]),etl3x=c(bootlist14[low],bootlist14[up]),rl3x=c(bootlist15[low],bootlist15[up]),ql3x=c(bootlist16[low],bootlist16[up]),
        
        skewx=c(bootlist17[low],bootlist17[up]),etskewx=c(bootlist18[low],bootlist18[up]),rskewx=c(bootlist19[low],bootlist19[up]),qskewx=c(bootlist20[low],bootlist20[up]),
        
        l4x=c(bootlist21[low],bootlist21[up]),etl4x=c(bootlist22[low],bootlist22[up]),rl4x=c(bootlist23[low],bootlist23[up]),ql4x=c(bootlist24[low],bootlist24[up]),
        
        kurtx=c(bootlist25[low],bootlist25[up]),etkurtx=c(bootlist26[low],bootlist26[up]),rkurtx=c(bootlist27[low],bootlist27[up]),qkurtx=c(bootlist28[low],bootlist28[up]),
        
        meany=c(bootlist29[low],bootlist29[up]),etmy=c(bootlist30[low],bootlist30[up]),rmy=c(bootlist31[low],bootlist31[up]),qmy=c(bootlist32[low],bootlist32[up]),
        l2y=c(bootlist33[low],bootlist33[up]),etl2y=c(bootlist34[low],bootlist34[up]),rl2y=c(bootlist35[low],bootlist35[up]),ql2y=c(bootlist36[low],bootlist36[up]),
        sdy=c(bootlist37[low],bootlist37[up]),etsdy=c(bootlist38[low],bootlist38[up]),
        rsdy=c(bootlist39[low],bootlist39[up]),qsdy=c(bootlist40[low],bootlist40[up]),
        l3y=c(bootlist41[low],bootlist41[up]),etl3y=c(bootlist42[low],bootlist42[up]),rl3y=c(bootlist43[low],bootlist43[up]),ql3y=c(bootlist44[low],bootlist44[up]),
        
        skewy=c(bootlist45[low],bootlist45[up]),etskewy=c(bootlist46[low],bootlist46[up]),rskewy=c(bootlist47[low],bootlist47[up]),qskewy=c(bootlist48[low],bootlist48[up]),
        
        l4y=c(bootlist49[low],bootlist49[up]),etl4y=c(bootlist50[low],bootlist50[up]),rl4y=c(bootlist51[low],bootlist51[up]),ql4y=c(bootlist52[low],bootlist52[up]),
        
        kurty=c(bootlist53[low],bootlist53[up]),etkurty=c(bootlist54[low],bootlist54[up]),rkurty=c(bootlist55[low],bootlist55[up]),qkurty=c(bootlist56[low],bootlist56[up])
  )
  
  se<-c(meanx=sd(bootlist1),etmx=sd(bootlist2),rmx=sd(bootlist3),qmx=sd(bootlist4),
        l2x=sd(bootlist5),etl2x=sd(bootlist6),rl2x=sd(bootlist7),ql2x=sd(bootlist8),
        sdx=sd(bootlist9),
        etsdx=sd(bootlist10),
        rsdx=sd(bootlist11),qsdx=sd(bootlist12),
        l3x=sd(bootlist13),etl3x=sd(bootlist14),rl3x=sd(bootlist15),ql3x=sd(bootlist16),
        
        skewx=sd(bootlist17),etskewx=sd(bootlist18),rskewx=sd(bootlist19),qskewx=sd(bootlist20),
        
        l4x=sd(bootlist21),etl4x=sd(bootlist22),rl4x=sd(bootlist23),ql4x=sd(bootlist24),
        
        kurtx=sd(bootlist25),etkurtx=sd(bootlist26),rkurtx=sd(bootlist27),qkurtx=sd(bootlist28),
        
        meany=sd(bootlist29),etmy=sd(bootlist30),rmy=sd(bootlist31),qmy=sd(bootlist32),
        l2y=sd(bootlist33),etl2y=sd(bootlist34),rl2y=sd(bootlist35),ql2y=sd(bootlist36),
        sdy=sd(bootlist37),etsdy=sd(bootlist38),
        rsdy=sd(bootlist39),qsdy=sd(bootlist40),
        l3y=sd(bootlist41),etl3y=sd(bootlist42),rl3y=sd(bootlist43),ql3y=sd(bootlist44),
        
        skewy=sd(bootlist45),etskewy=sd(bootlist46),rskewy=sd(bootlist47),qskewy=sd(bootlist48),
        
        l4y=sd(bootlist49),etl4y=sd(bootlist50),rl4y=sd(bootlist51),ql4y=sd(bootlist52),
        
        kurty=sd(bootlist53),etkurty=sd(bootlist54),rkurty=sd(bootlist55),qkurty=sd(bootlist56)
  )
  
  all<-list(p_value_diff=p_value_diff,ci_diff=ci_diff,ci=ci,se=se,estimate=estimate1)
  if((rqfmx[8])/((rqscalex[8])^(2))<4.7 & standist=="exponential"){
    print("The quantile kurtosis is lower than 4.7, it might be better to use the Rayleigh distribution as the standard distribution.")
  }
  return(all)
}

ebh2parallel<-function (x,y,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray"),alpha=0.05,nboot=100){
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
    mmm<-function(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      etm1<-etm(sortedx,interval=interval,fast=fast,batch=batch)
      if(etm1[2]==Inf){
        return(print("ETM is infinity, due to the double precision floating point limits. Usually, the solution is transforming your original data."))
      }
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
      output1<-c(mean=mean(sortedx),etm=etm1[2],rm=rm1,qm=qm1)
      return(output1)
    }
    rqmean<-function (x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = FALSE,alpha = 0.05,nboot = 1000){
      if(cise){
        return (mmmci(x, interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm,alpha = alpha,nboot = nboot))
      } 
      else {return (mmm(x,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm))}
    }
    
    rqscale<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.3659,drm=0.7930,dlqm=0.8218,dqm=0.7825,sd=FALSE){
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
      lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      if(sd){
        lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drm,dqm=drm)
        msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drm,dqm=drm)
        allmo<-c(l2=lmo[1],etl2=lmo[2],rl2=lmo[3],ql2=lmo[4],
                 var=mo[1],etvar=mo[2],rvar=mo[3],qvar=mo[4])
        allsd<-c(l2sd=lmsd[1],etl2sd=lmsd[2],rl2sd=lmsd[3],ql2sd=lmsd[4],
                 varsd=msd[1],etvarsd=msd[2],rvarsd=msd[3],qvarsd=msd[4])
        all<-c(allmo,allsd)
      }else{
        all<-c(l2=lmo[1],etl2=lmo[2],rl2=lmo[3],ql2=lmo[4],
               var=mo[1],etvar=mo[2],rvar=mo[3],qvar=mo[4])
      }
      return(all)
    }
    
    rqsd<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,drm=0.7930,dqm=0.7825){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthn<-length(sortedx)
      if (boot){
        subtract<-t(replicate(times , sort(sample(sortedx, size = 2))))
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
        dp2m<-(dp^2)/2
      }
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      all<-c(sd=sqrt(mo[1]),etsd=sqrt(mo[2]),rsd=sqrt(mo[3]),qsd=sqrt(mo[4]))
      return(all)
    }
    
    rqtm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.1808,drm=1.7492,dlqm=1.1753,dqm=0.5715,sd=FALSE,drsd=0.7930,dqsd=0.7825){
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
      lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      if(sd){
        lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        allmo<-c(l3=lmo[1],etl3=lmo[2],rl3=lmo[3],ql3=lmo[4],
                 tm=mo[1],ettm=mo[2],rtm=mo[3],qtm=mo[4])
        allsd<-c(l3sd=lmsd[1],etl3sd=lmsd[2],rl3sd=lmsd[3],ql3sd=lmsd[4],
                 tmsd=msd[1],ettmsd=msd[2],rtmsd=msd[3],qtmsd=msd[4])
        all<-c(allmo,allsd)
      }else{
        all<-c(l3=lmo[1],etl3=lmo[2],rl3=lmo[3],ql3=lmo[4],
               tm=mo[1],ettm=mo[2],rtm=mo[3],qtm=mo[4])
      }
      return(all)
    }
    
    rqfm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=-0.3542,drm=3.4560,dlqm=NaN,dqm=0.1246,sd=FALSE,drsd=0.7930,dqsd=0.7825){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthn<-length(sortedx)
      
      getm<-function(vector){ 
        resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                      4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                      vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                      6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                      vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                                   6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
        return(resd)
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
      
      lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      if(sd){
        lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        allmo<-c(l4=lmo[1],etl4=lmo[2],rl4=lmo[3],ql4=lmo[4],
                 fm=mo[1],etfm=mo[2],rfm=mo[3],qfm=mo[4])
        allsd<-c(l4sd=lmsd[1],etl4sd=lmsd[2],rl4sd=lmsd[3],ql4sd=lmsd[4],
                 fmsd=msd[1],etfmsd=msd[2],rfmsd=msd[3],qfmsd=msd[4])
        all<-c(allmo,allsd)
      }else{
        all<-c(l4=lmo[1],etl4=lmo[2],rl4=lmo[3],ql4=lmo[4],
               fm=mo[1],etfm=mo[2],rfm=mo[3],qfm=mo[4])
      }
      return(all)
    }
    mmmx<-mmm(x=datax[i,],interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
    rqscalex<-rqscale(x=datax[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=FALSE)
    rqtmx<-rqtm(x=datax[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
    rqfmx<-rqfm(x=datax[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
    mmmy<-mmm(x=datay[i,],interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
    rqscaley<-rqscale(x=datay[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=FALSE)
    rqtmy<-rqtm(x=datay[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
    rqfmy<-rqfm(x=datay[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
    
    estimate1<-c(c(meanx=mmmx[1],etmx=mmmx[2],rmx=mmmx[3],qmx=mmmx[4]),
                 c(l2x=rqscalex[1],etl2x=rqscalex[2],rl2x=rqscalex[3],ql2x=rqscalex[4],sdx=sqrt(rqscalex[5]),etsdx=sqrt(rqscalex[6]),
                   rsdx=sqrt(rqscalex[7]),qsdx=sqrt(rqscalex[8])),
                 c(l3x=rqtmx[1]/rqscalex[1],etl3x=rqtmx[2]/rqscalex[2],rl3x=rqtmx[3]/rqscalex[3],ql3x=rqtmx[4]/rqscalex[4],
                   skewx=(rqtmx[5])/((rqscalex[5])^(3/2)),etskewx=(rqtmx[6])/((rqscalex[6])^(3/2)),
                   rskewx=(rqtmx[7])/((rqscalex[7])^(3/2)),qskewx=(rqtmx[8])/((rqscalex[8])^(3/2))),
                 c(l4x=rqfmx[1]/rqscalex[1],etl4x=rqfmx[2]/rqscalex[2],rl4x=rqfmx[3]/rqscalex[3],ql4x=rqfmx[4]/rqscalex[4],kurtx=(rqfmx[5])/((rqscalex[5])^(2)),etkurtx=(rqfmx[6])/((rqscalex[6])^(2)),rkurtx=(rqfmx[7])/((rqscalex[7])^(2)),qkurtx=(rqfmx[8])/((rqscalex[8])^(2))),
                 c(meany=mmmy[1],etmy=mmmy[2],rmy=mmmy[3],qmy=mmmy[4]),
                 c(l2y=rqscaley[1],etl2y=rqscaley[2],rl2y=rqscaley[3],ql2y=rqscaley[4],sdy=sqrt(rqscaley[5]),etsdy=sqrt(rqscaley[6]),
                   rsdy=sqrt(rqscaley[7]),qsdy=sqrt(rqscaley[8])),
                 c(l3y=rqtmy[1]/rqscaley[1],etl3y=rqtmy[2]/rqscaley[2],rl3y=rqtmy[3]/rqscaley[3],ql3y=rqtmy[4]/rqscaley[4],
                   skewy=(rqtmy[5])/((rqscaley[5])^(3/2)),etskewy=(rqtmy[6])/((rqscaley[6])^(3/2)),
                   rskewy=(rqtmy[7])/((rqscaley[7])^(3/2)),qskewy=(rqtmy[8])/((rqscaley[8])^(3/2))),
                 c(l4y=rqfmy[1]/rqscaley[1],etl4y=rqfmy[2]/rqscaley[2],rl4y=rqfmy[3]/rqscaley[3],ql4y=rqfmy[4]/rqscaley[4],kurty=(rqfmy[5])/((rqscaley[5])^(2)),etkurty=(rqfmy[6])/((rqscaley[6])^(2)),rkurty=(rqfmy[7])/((rqscaley[7])^(2)),qkurty=(rqfmy[8])/((rqscaley[8])^(2)))
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
  bootlist45<-(as.matrix(estimate0[,45]))
  bootlist46<-(as.matrix(estimate0[,46]))
  bootlist47<-(as.matrix(estimate0[,47]))
  bootlist48<-(as.matrix(estimate0[,48]))
  bootlist49<-(as.matrix(estimate0[,49]))
  bootlist50<-(as.matrix(estimate0[,50]))
  bootlist51<-(as.matrix(estimate0[,51]))
  bootlist52<-(as.matrix(estimate0[,52]))
  bootlist53<-(as.matrix(estimate0[,53]))
  bootlist54<-(as.matrix(estimate0[,54]))
  bootlist55<-(as.matrix(estimate0[,55]))
  bootlist56<-(as.matrix(estimate0[,56]))
  
  bootlist1a<-sort(bootlist1-bootlist29)
  bootlist2a<-sort(bootlist2-bootlist30)
  bootlist3a<-sort(bootlist3-bootlist31)
  bootlist4a<-sort(bootlist4-bootlist32)
  bootlist5a<-sort(bootlist5-bootlist33)
  bootlist6a<-sort(bootlist6-bootlist34)
  bootlist7a<-sort(bootlist7-bootlist35)
  bootlist8a<-sort(bootlist8-bootlist36)
  bootlist9a<-sort(bootlist9-bootlist37)
  bootlist10a<-sort(bootlist10-bootlist38)
  bootlist11a<-sort(bootlist11-bootlist39)
  bootlist12a<-sort(bootlist12-bootlist40)
  bootlist13a<-sort(bootlist13-bootlist41)
  bootlist14a<-sort(bootlist14-bootlist42)
  bootlist15a<-sort(bootlist15-bootlist43)
  bootlist16a<-sort(bootlist16-bootlist44)
  bootlist17a<-sort(bootlist17-bootlist45)
  bootlist18a<-sort(bootlist18-bootlist46)
  bootlist19a<-sort(bootlist19-bootlist47)
  bootlist20a<-sort(bootlist20-bootlist48)
  bootlist21a<-sort(bootlist21-bootlist49)
  bootlist22a<-sort(bootlist22-bootlist50)
  bootlist23a<-sort(bootlist23-bootlist51)
  bootlist24a<-sort(bootlist24-bootlist52)
  bootlist25a<-sort(bootlist25-bootlist53)
  bootlist26a<-sort(bootlist26-bootlist54)
  bootlist27a<-sort(bootlist27-bootlist55)
  bootlist28a<-sort(bootlist28-bootlist56)
  
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
  bootlist45<-sort(bootlist45)
  bootlist46<-sort(bootlist46)
  bootlist47<-sort(bootlist47)
  bootlist48<-sort(bootlist48)
  bootlist49<-sort(bootlist49)
  bootlist50<-sort(bootlist50)
  bootlist51<-sort(bootlist51)
  bootlist52<-sort(bootlist52)
  bootlist53<-sort(bootlist53)
  bootlist54<-sort(bootlist54)
  bootlist55<-sort(bootlist55)
  bootlist56<-sort(bootlist56)
  
  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  
  mmmx<-mmm(x=sortedx,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  rqscalex<-rqscale(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=FALSE)
  rqtmx<-rqtm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
  rqfmx<-rqfm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
  mmmy<-mmm(x=sortedy,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  rqscaley<-rqscale(x=sortedy,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=FALSE)
  rqtmy<-rqtm(x=sortedy,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
  rqfmy<-rqfm(x=sortedy,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=FALSE,drsd=drmscale,dqsd=dqmscale)
  
  estimate<-c(c(meanx=mmmx[1],etmx=mmmx[2],rmx=mmmx[3],qmx=mmmx[4]),
               c(l2x=rqscalex[1],etl2x=rqscalex[2],rl2x=rqscalex[3],ql2x=rqscalex[4],sdx=sqrt(rqscalex[5]),etsdx=sqrt(rqscalex[6]),
                 rsdx=sqrt(rqscalex[7]),qsdx=sqrt(rqscalex[8])),
               c(l3x=rqtmx[1]/rqscalex[1],etl3x=rqtmx[2]/rqscalex[2],rl3x=rqtmx[3]/rqscalex[3],ql3x=rqtmx[4]/rqscalex[4],
                 skewx=(rqtmx[5])/((rqscalex[5])^(3/2)),etskewx=(rqtmx[6])/((rqscalex[6])^(3/2)),
                 rskewx=(rqtmx[7])/((rqscalex[7])^(3/2)),qskewx=(rqtmx[8])/((rqscalex[8])^(3/2))),
               c(l4x=rqfmx[1]/rqscalex[1],etl4x=rqfmx[2]/rqscalex[2],rl4x=rqfmx[3]/rqscalex[3],ql4x=rqfmx[4]/rqscalex[4],kurtx=(rqfmx[5])/((rqscalex[5])^(2)),etkurtx=(rqfmx[6])/((rqscalex[6])^(2)),rkurtx=(rqfmx[7])/((rqscalex[7])^(2)),qkurtx=(rqfmx[8])/((rqscalex[8])^(2))),
               c(meany=mmmy[1],etmy=mmmy[2],rmy=mmmy[3],qmy=mmmy[4]),
               c(l2y=rqscaley[1],etl2y=rqscaley[2],rl2y=rqscaley[3],ql2y=rqscaley[4],sdy=sqrt(rqscaley[5]),etsdy=sqrt(rqscaley[6]),
                 rsdy=sqrt(rqscaley[7]),qsdy=sqrt(rqscaley[8])),
               c(l3y=rqtmy[1]/rqscaley[1],etl3y=rqtmy[2]/rqscaley[2],rl3y=rqtmy[3]/rqscaley[3],ql3y=rqtmy[4]/rqscaley[4],
                 skewy=(rqtmy[5])/((rqscaley[5])^(3/2)),etskewy=(rqtmy[6])/((rqscaley[6])^(3/2)),
                 rskewy=(rqtmy[7])/((rqscaley[7])^(3/2)),qskewy=(rqtmy[8])/((rqscaley[8])^(3/2))),
               c(l4y=rqfmy[1]/rqscaley[1],etl4y=rqfmy[2]/rqscaley[2],rl4y=rqfmy[3]/rqscaley[3],ql4y=rqfmy[4]/rqscaley[4],kurty=(rqfmy[5])/((rqscaley[5])^(2)),etkurty=(rqfmy[6])/((rqscaley[6])^(2)),rkurty=(rqfmy[7])/((rqscaley[7])^(2)),qkurty=(rqfmy[8])/((rqscaley[8])^(2)))
  )
  
  
  bootlist1a<-bootlist1a-(estimate[1]-estimate[29])
  bootlist2a<-bootlist2a-(estimate[2]-estimate[30])
  bootlist3a<-bootlist3a-(estimate[3]-estimate[31])
  bootlist4a<-bootlist4a-(estimate[4]-estimate[32])
  bootlist5a<-bootlist5a-(estimate[5]-estimate[33])
  bootlist6a<-bootlist6a-(estimate[6]-estimate[34])
  bootlist7a<-bootlist7a-(estimate[7]-estimate[35])
  bootlist8a<-bootlist8a-(estimate[8]-estimate[36])
  bootlist9a<-bootlist9a-(estimate[9]-estimate[37])
  bootlist10a<-bootlist10a-(estimate[10]-estimate[38])
  bootlist11a<-bootlist11a-(estimate[11]-estimate[39])
  bootlist12a<-bootlist12a-(estimate[12]-estimate[40])
  bootlist13a<-bootlist13a-(estimate[13]-estimate[41])
  bootlist14a<-bootlist14a-(estimate[14]-estimate[42])
  bootlist15a<-bootlist15a-(estimate[15]-estimate[43])
  bootlist16a<-bootlist16a-(estimate[16]-estimate[44])
  bootlist17a<-bootlist17a-(estimate[17]-estimate[45])
  bootlist18a<-bootlist18a-(estimate[18]-estimate[46])
  bootlist19a<-bootlist19a-(estimate[19]-estimate[47])
  bootlist20a<-bootlist20a-(estimate[20]-estimate[48])
  bootlist21a<-bootlist21a-(estimate[21]-estimate[49])
  bootlist22a<-bootlist22a-(estimate[22]-estimate[50])
  bootlist23a<-bootlist23a-(estimate[23]-estimate[51])
  bootlist24a<-bootlist24a-(estimate[24]-estimate[52])
  bootlist25a<-bootlist25a-(estimate[25]-estimate[53])
  bootlist26a<-bootlist26a-(estimate[26]-estimate[54])
  bootlist27a<-bootlist27a-(estimate[27]-estimate[55])
  bootlist28a<-bootlist28a-(estimate[28]-estimate[56])
  
  cidiff<-function(bootlist,low,up,estimate,n){
    (estimate[n]-estimate[n+28])+c(bootlist[low],bootlist[up])
  }
  pb2<-function(bootlist,estimate,n){
    null_value<-(estimate[n]-estimate[n+28])
    p <- mean(bootlist < null_value) + 0.5 * mean(bootlist == null_value)
    p <- 2 * min(c(p, 1 - p))
    p
  }
  ci_diff<-c(c(mean=cidiff(bootlist1a,low=low,up=up,estimate=estimate,n=1),etm=cidiff(bootlist2a,low=low,up=up,estimate=estimate,n=2),rm=cidiff(bootlist3a,low=low,up=up,estimate=estimate,n=3),qm=cidiff(bootlist4a,low=low,up=up,estimate=estimate,n=4)),
             c(l2=cidiff(bootlist5a,low=low,up=up,estimate=estimate,n=5),etl2=cidiff(bootlist6a,low=low,up=up,estimate=estimate,n=6),rl2=cidiff(bootlist7a,low=low,up=up,estimate=estimate,n=7),ql2=cidiff(bootlist8a,low=low,up=up,estimate=estimate,n=8),sd=cidiff(bootlist9a,low=low,up=up,estimate=estimate,n=9),etsd=cidiff(bootlist10a,low=low,up=up,estimate=estimate,n=10),
               rsd=cidiff(bootlist11a,low=low,up=up,estimate=estimate,n=11),qsd=cidiff(bootlist12a,low=low,up=up,estimate=estimate,n=12)),
             c(l3=cidiff(bootlist13a,low=low,up=up,estimate=estimate,n=13),etl3=cidiff(bootlist14a,low=low,up=up,estimate=estimate,n=14),rl3=cidiff(bootlist15a,low=low,up=up,estimate=estimate,n=15),ql3=cidiff(bootlist16a,low=low,up=up,estimate=estimate,n=16),
               skew=cidiff(bootlist17a,low=low,up=up,estimate=estimate,n=17),etskew=cidiff(bootlist18a,low=low,up=up,estimate=estimate,n=18),
               rskew=cidiff(bootlist19a,low=low,up=up,estimate=estimate,n=19),qskew=cidiff(bootlist20a,low=low,up=up,estimate=estimate,n=20)),
             c(l4=cidiff(bootlist21a,low=low,up=up,estimate=estimate,n=21),etl4=cidiff(bootlist22a,low=low,up=up,estimate=estimate,n=22),rl4=cidiff(bootlist23a,low=low,up=up,estimate=estimate,n=23),ql4=cidiff(bootlist24a,low=low,up=up,estimate=estimate,n=24),kurt=cidiff(bootlist25a,low=low,up=up,estimate=estimate,n=25),etkurt=cidiff(bootlist26a,low=low,up=up,estimate=estimate,n=26),rkurt=cidiff(bootlist27a,low=low,up=up,estimate=estimate,n=27),qkurt=cidiff(bootlist28a,low=low,up=up,estimate=estimate,n=28)))
  
  p_value_diff<-c(c(mean=pb2(bootlist1a,estimate=estimate,n=1),etm=pb2(bootlist2a,estimate=estimate,n=2),rm=pb2(bootlist3a,estimate=estimate,n=3),qm=pb2(bootlist4a,estimate=estimate,n=4)),
                  c(l2=pb2(bootlist5a,estimate=estimate,n=5),etl2=pb2(bootlist6a,estimate=estimate,n=6),rl2=pb2(bootlist7a,estimate=estimate,n=7),ql2=pb2(bootlist8a,estimate=estimate,n=8),sd=pb2(bootlist9a,estimate=estimate,n=9),etsd=pb2(bootlist10a,estimate=estimate,n=10),
                    rsd=pb2(bootlist11a,estimate=estimate,n=11),qsd=pb2(bootlist12a,estimate=estimate,n=12)),
                  c(l3=pb2(bootlist13a,estimate=estimate,n=13),etl3=pb2(bootlist14a,estimate=estimate,n=14),rl3=pb2(bootlist15a,estimate=estimate,n=15),ql3=pb2(bootlist16a,estimate=estimate,n=16),
                    skew=pb2(bootlist17a,estimate=estimate,n=17),etskew=pb2(bootlist18a,estimate=estimate,n=18),
                    rskew=pb2(bootlist19a,estimate=estimate,n=19),qskew=pb2(bootlist20a,estimate=estimate,n=20)),
                  c(l4=pb2(bootlist21a,estimate=estimate,n=21),etl4=pb2(bootlist22a,estimate=estimate,n=22),rl4=pb2(bootlist23a,estimate=estimate,n=23),ql4=pb2(bootlist24a,estimate=estimate,n=24),kurt=pb2(bootlist25a,estimate=estimate,n=25),etkurt=pb2(bootlist26a,estimate=estimate,n=26),rkurt=pb2(bootlist27a,estimate=estimate,n=27),qkurt=pb2(bootlist28a,estimate=estimate,n=28)))
  
  
  ci<-c(meanx=c(bootlist1[low],bootlist1[up]),etmx=c(bootlist2[low],bootlist2[up]),rmx=c(bootlist3[low],bootlist3[up]),qmx=c(bootlist4[low],bootlist4[up]),
        l2x=c(bootlist5[low],bootlist5[up]),etl2x=c(bootlist6[low],bootlist6[up]),rl2x=c(bootlist7[low],bootlist7[up]),ql2x=c(bootlist8[low],bootlist8[up]),
        sdx=c(bootlist9[low],bootlist9[up]),
        etsdx=c(bootlist10[low],bootlist10[up]),
        rsdx=c(bootlist11[low],bootlist11[up]),qsdx=c(bootlist12[low],bootlist12[up]),
        l3x=c(bootlist13[low],bootlist13[up]),etl3x=c(bootlist14[low],bootlist14[up]),rl3x=c(bootlist15[low],bootlist15[up]),ql3x=c(bootlist16[low],bootlist16[up]),
        
        skewx=c(bootlist17[low],bootlist17[up]),etskewx=c(bootlist18[low],bootlist18[up]),rskewx=c(bootlist19[low],bootlist19[up]),qskewx=c(bootlist20[low],bootlist20[up]),
        
        l4x=c(bootlist21[low],bootlist21[up]),etl4x=c(bootlist22[low],bootlist22[up]),rl4x=c(bootlist23[low],bootlist23[up]),ql4x=c(bootlist24[low],bootlist24[up]),
        
        kurtx=c(bootlist25[low],bootlist25[up]),etkurtx=c(bootlist26[low],bootlist26[up]),rkurtx=c(bootlist27[low],bootlist27[up]),qkurtx=c(bootlist28[low],bootlist28[up]),
        
        meany=c(bootlist29[low],bootlist29[up]),etmy=c(bootlist30[low],bootlist30[up]),rmy=c(bootlist31[low],bootlist31[up]),qmy=c(bootlist32[low],bootlist32[up]),
        l2y=c(bootlist33[low],bootlist33[up]),etl2y=c(bootlist34[low],bootlist34[up]),rl2y=c(bootlist35[low],bootlist35[up]),ql2y=c(bootlist36[low],bootlist36[up]),
        sdy=c(bootlist37[low],bootlist37[up]),etsdy=c(bootlist38[low],bootlist38[up]),
        rsdy=c(bootlist39[low],bootlist39[up]),qsdy=c(bootlist40[low],bootlist40[up]),
        l3y=c(bootlist41[low],bootlist41[up]),etl3y=c(bootlist42[low],bootlist42[up]),rl3y=c(bootlist43[low],bootlist43[up]),ql3y=c(bootlist44[low],bootlist44[up]),
        
        skewy=c(bootlist45[low],bootlist45[up]),etskewy=c(bootlist46[low],bootlist46[up]),rskewy=c(bootlist47[low],bootlist47[up]),qskewy=c(bootlist48[low],bootlist48[up]),
        
        l4y=c(bootlist49[low],bootlist49[up]),etl4y=c(bootlist50[low],bootlist50[up]),rl4y=c(bootlist51[low],bootlist51[up]),ql4y=c(bootlist52[low],bootlist52[up]),
        
        kurty=c(bootlist53[low],bootlist53[up]),etkurty=c(bootlist54[low],bootlist54[up]),rkurty=c(bootlist55[low],bootlist55[up]),qkurty=c(bootlist56[low],bootlist56[up])
  )
  
  se<-c(meanx=sd(bootlist1),etmx=sd(bootlist2),rmx=sd(bootlist3),qmx=sd(bootlist4),
        l2x=sd(bootlist5),etl2x=sd(bootlist6),rl2x=sd(bootlist7),ql2x=sd(bootlist8),
        sdx=sd(bootlist9),
        etsdx=sd(bootlist10),
        rsdx=sd(bootlist11),qsdx=sd(bootlist12),
        l3x=sd(bootlist13),etl3x=sd(bootlist14),rl3x=sd(bootlist15),ql3x=sd(bootlist16),
        
        skewx=sd(bootlist17),etskewx=sd(bootlist18),rskewx=sd(bootlist19),qskewx=sd(bootlist20),
        
        l4x=sd(bootlist21),etl4x=sd(bootlist22),rl4x=sd(bootlist23),ql4x=sd(bootlist24),
        
        kurtx=sd(bootlist25),etkurtx=sd(bootlist26),rkurtx=sd(bootlist27),qkurtx=sd(bootlist28),
        
        meany=sd(bootlist29),etmy=sd(bootlist30),rmy=sd(bootlist31),qmy=sd(bootlist32),
        l2y=sd(bootlist33),etl2y=sd(bootlist34),rl2y=sd(bootlist35),ql2y=sd(bootlist36),
        sdy=sd(bootlist37),etsdy=sd(bootlist38),
        rsdy=sd(bootlist39),qsdy=sd(bootlist40),
        l3y=sd(bootlist41),etl3y=sd(bootlist42),rl3y=sd(bootlist43),ql3y=sd(bootlist44),
        
        skewy=sd(bootlist45),etskewy=sd(bootlist46),rskewy=sd(bootlist47),qskewy=sd(bootlist48),
        
        l4y=sd(bootlist49),etl4y=sd(bootlist50),rl4y=sd(bootlist51),ql4y=sd(bootlist52),
        
        kurty=sd(bootlist53),etkurty=sd(bootlist54),rkurty=sd(bootlist55),qkurty=sd(bootlist56)
  )
  
  all<-list(p_value_diff=p_value_diff,ci_diff=ci_diff,ci=ci,se=se,estimate=estimate)
  if((rqfmx[8])/((rqscalex[8])^(2))<4.7 & standist=="exponential"){
    print("The quantile kurtosis is lower than 4.7, it might be better to use the Rayleigh distribution as the standard distribution.")
  }
  return(all)
}

NRSs<-function(x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray"),sd=FALSE,cise = FALSE,parallel=TRUE,alpha = 0.05,nboot = 100,null_mean=1,null_sd=1,null_skew=2,null_kurt=9,null_l2=0.5,null_l3=1/3,null_l4=1/6){
  if (times%%9!=0){
    return ("Please set times as a multiple of 9.")
  }
  if(cise & parallel){
    return (NRSsciparallel(x, interval=interval,fast=fast,batch=batch,boot=boot,times =times ,standist=standist,alpha=alpha,nboot=nboot,null_mean=null_mean,null_sd=null_sd,null_skew=null_skew,null_kurt=null_kurt,null_l2=null_l2,null_l3=null_l3,null_l4=null_l4))
  } else if(cise){return (NRSsci(x, interval=interval,fast=fast,batch=batch,boot=boot,times =times ,standist=standist,alpha=alpha,nboot=nboot,null_mean=null_mean,null_sd=null_sd,null_skew=null_skew,null_kurt=null_kurt,null_l2=null_l2,null_l3=null_l3,null_l4=null_l4))
  }
  else{return (NRSssimple(x, interval=interval,fast=fast,batch=batch,boot=boot,times =times ,standist=standist,sd=sd))
}}

htest<-function(x,y,boottype=c("empirial","percentile"),interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray"),alpha=0.05,nboot=100){
  if (boottype=="empirial"){
    return(ebh2parallel(x,y,interval=interval,fast=fast,batch=batch,boot=boot,times =times,standist=standist,alpha=alpha,nboot=nboot))
  } 
  else{return (pbh2parallel(x,y,interval=interval,fast=fast,batch=batch,boot=boot,times =times,standist=standist,alpha=alpha,nboot=nboot))
}}
effectsizeNRSssimple<-function(x,y,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray")){
  sortedx<-sort(x,decreasing = FALSE,method ="radix")
  lengthx<-length(sortedx)
  sortedy<-sort(y,decreasing = FALSE,method ="radix")
  lengthy<-length(sortedy)
  estimatex<-NRSssimple(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,standist=standist,sd=TRUE)
  estimatey<-NRSssimple(x=sortedy,interval=interval,fast=fast,batch=batch,boot=boot,times =times,standist=standist,sd=TRUE)
  firsteffectsize<-c(mean=(estimatex$first[1]-estimatey$first[1])/(((((estimatex$first[1])^2)+((estimatey$first[1])^2))*0.5)^0.5),etm=(estimatex$first[2]-estimatey$first[2])/(((((estimatex$first[2])^2)+((estimatey$first[2])^2))*0.5)^0.5),rm=(estimatex$first[3]-estimatey$first[3])/(((((estimatex$first[3])^2)+((estimatey$first[3])^2))*0.5)^0.5),qm=(estimatex$first[4]-estimatey$first[4])/(((((estimatex$first[4])^2)+((estimatey$first[4])^2))*0.5)^0.5))
  
  secondeffectsize<-c(l2=(estimatex$second[1]-estimatey$second[1])/(((((estimatex$second[1])^2)+((estimatey$second[1])^2))*0.5)^0.5),etl2=(estimatex$second[2]-estimatey$second[2])/(((((estimatex$second[2])^2)+((estimatey$second[2])^2))*0.5)^0.5),rl2=(estimatex$second[3]-estimatey$second[3])/(((((estimatex$second[3])^2)+((estimatey$second[3])^2))*0.5)^0.5),ql2=(estimatex$second[4]-estimatey$second[4])/(((((estimatex$second[4])^2)+((estimatey$second[4])^2))*0.5)^0.5),
                      sd=(estimatex$second[5]-estimatey$second[5])/(((((estimatex$second[5])^2)+((estimatey$second[5])^2))*0.5)^0.5),etsd=(estimatex$second[6]-estimatey$second[6])/(((((estimatex$second[6])^2)+((estimatey$second[6])^2))*0.5)^0.5),rsd=(estimatex$second[7]-estimatey$second[7])/(((((estimatex$second[7])^2)+((estimatey$second[7])^2))*0.5)^0.5),qsd=(estimatex$second[8]-estimatey$second[8])/(((((estimatex$second[8])^2)+((estimatey$second[8])^2))*0.5)^0.5)
  )
  thirdeffectsize<-c(l3=(estimatex$third[1]-estimatey$third[1])/(((((estimatex$third[1])^2)+((estimatey$third[1])^2))*0.5)^0.5),etl3=(estimatex$third[2]-estimatey$third[2])/(((((estimatex$third[2])^2)+((estimatey$third[2])^2))*0.5)^0.5),rl3=(estimatex$third[3]-estimatey$third[3])/(((((estimatex$third[3])^2)+((estimatey$third[3])^2))*0.5)^0.5),ql3=(estimatex$third[4]-estimatey$third[4])/(((((estimatex$third[4])^2)+((estimatey$third[4])^2))*0.5)^0.5),
                     skew=(estimatex$third[5]-estimatey$third[5])/(((((estimatex$third[5])^2)+((estimatey$third[5])^2))*0.5)^0.5),etskew=(estimatex$third[6]-estimatey$third[6])/(((((estimatex$third[6])^2)+((estimatey$third[6])^2))*0.5)^0.5),rskew=(estimatex$third[7]-estimatey$third[7])/(((((estimatex$third[7])^2)+((estimatey$third[7])^2))*0.5)^0.5),qskew=(estimatex$third[8]-estimatey$third[8])/(((((estimatex$third[8])^2)+((estimatey$third[8])^2))*0.5)^0.5)
  )
  fourtheffectsize<-c(l4=(estimatex$fourth[1]-estimatey$fourth[1])/(((((estimatex$fourth[1])^2)+((estimatey$fourth[1])^2))*0.5)^0.5),etl4=(estimatex$fourth[2]-estimatey$fourth[2])/(((((estimatex$fourth[2])^2)+((estimatey$fourth[2])^2))*0.5)^0.5),rl4=(estimatex$fourth[3]-estimatey$fourth[3])/(((((estimatex$fourth[3])^2)+((estimatey$fourth[3])^2))*0.5)^0.5),ql4=(estimatex$fourth[4]-estimatey$fourth[4])/(((((estimatex$fourth[4])^2)+((estimatey$fourth[4])^2))*0.5)^0.5),
                      kurt=(estimatex$fourth[5]-estimatey$fourth[5])/(((((estimatex$fourth[5])^2)+((estimatey$fourth[5])^2))*0.5)^0.5),etkurt=(estimatex$fourth[6]-estimatey$fourth[6])/(((((estimatex$fourth[6])^2)+((estimatey$fourth[6])^2))*0.5)^0.5),rkurt=(estimatex$fourth[7]-estimatey$fourth[7])/(((((estimatex$fourth[7])^2)+((estimatey$fourth[7])^2))*0.5)^0.5),qkurt=(estimatex$fourth[8]-estimatey$fourth[8])/(((((estimatex$fourth[8])^2)+((estimatey$fourth[8])^2))*0.5)^0.5)
  )
  
  all<-list(firsteffectsize=firsteffectsize,secondeffectsize=secondeffectsize,thirdeffectsize=thirdeffectsize,fourtheffectsize=fourtheffectsize,estimatex=estimatex,estimatey=estimatey)
  if(estimatex$fourth[5]<4.7 & standist=="exponential"){
    print("The quantile kurtosis is lower than 4.7, it might be better to use the Rayleigh distribution as the standard distribution.")
  }
  
  return(all)
}

esbootparallel<-function (x,y,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray"),alpha=0.05,nboot=100){
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
    mmm<-function(x,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      etm1<-etm(sortedx,interval=interval,fast=fast,batch=batch)
      if(etm1[2]==Inf){
        return(print("ETM is infinity, due to the double precision floating point limits. Usually, the solution is transforming your original data."))
      }
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
      output1<-c(mean=mean(sortedx),etm=etm1[2],rm=rm1,qm=qm1)
      return(output1)
    }

    rqscale<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.3659,drm=0.7930,dlqm=0.8218,dqm=0.7825,sd=FALSE){
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
      lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      if(sd){
        lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drm,dqm=drm)
        msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drm,dqm=drm)
        allmo<-c(l2=lmo[1],etl2=lmo[2],rl2=lmo[3],ql2=lmo[4],
                 var=mo[1],etvar=mo[2],rvar=mo[3],qvar=mo[4])
        allsd<-c(l2sd=lmsd[1],etl2sd=lmsd[2],rl2sd=lmsd[3],ql2sd=lmsd[4],
                 varsd=msd[1],etvarsd=msd[2],rvarsd=msd[3],qvarsd=msd[4])
        all<-c(allmo,allsd)
      }else{
        all<-c(l2=lmo[1],etl2=lmo[2],rl2=lmo[3],ql2=lmo[4],
               var=mo[1],etvar=mo[2],rvar=mo[3],qvar=mo[4])
      }
      return(all)
    }
    
    rqsd<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,drm=0.7930,dqm=0.7825){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthn<-length(sortedx)
      if (boot){
        subtract<-t(replicate(times , sort(sample(sortedx, size = 2))))
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
        dp2m<-(dp^2)/2
      }
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      all<-c(sd=sqrt(mo[1]),etsd=sqrt(mo[2]),rsd=sqrt(mo[3]),qsd=sqrt(mo[4]))
      return(all)
    }
    
    rqtm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=0.1808,drm=1.7492,dlqm=1.1753,dqm=0.5715,sd=FALSE,drsd=0.7930,dqsd=0.7825){
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
      lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      if(sd){
        lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        allmo<-c(l3=lmo[1],etl3=lmo[2],rl3=lmo[3],ql3=lmo[4],
                 tm=mo[1],ettm=mo[2],rtm=mo[3],qtm=mo[4])
        allsd<-c(l3sd=lmsd[1],etl3sd=lmsd[2],rl3sd=lmsd[3],ql3sd=lmsd[4],
                 tmsd=msd[1],ettmsd=msd[2],rtmsd=msd[3],qtmsd=msd[4])
        all<-c(allmo,allsd)
      }else{
        all<-c(l3=lmo[1],etl3=lmo[2],rl3=lmo[3],ql3=lmo[4],
               tm=mo[1],ettm=mo[2],rtm=mo[3],qtm=mo[4])
      }
      return(all)
    }
    
    rqfm<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,dlrm=-0.3542,drm=3.4560,dlqm=NaN,dqm=0.1246,sd=FALSE,drsd=0.7930,dqsd=0.7825){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthn<-length(sortedx)
      
      getm<-function(vector){ 
        resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                      4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                      vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                      6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                      vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                                   6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
        return(resd)
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
      
      lmo<-mmm(x=dp2lm,interval=interval,fast=fast,batch=batch,drm=dlrm,dqm=dlqm)
      mo<-mmm(x=dp2m,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      if(sd){
        lmsd<-rqsd(x=dp2lm,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        msd<-rqsd(x=dp2m,interval=interval,fast=fast,batch=batch,boot=boot,times =times,drm=drsd,dqm=dqsd)
        allmo<-c(l4=lmo[1],etl4=lmo[2],rl4=lmo[3],ql4=lmo[4],
                 fm=mo[1],etfm=mo[2],rfm=mo[3],qfm=mo[4])
        allsd<-c(l4sd=lmsd[1],etl4sd=lmsd[2],rl4sd=lmsd[3],ql4sd=lmsd[4],
                 fmsd=msd[1],etfmsd=msd[2],rfmsd=msd[3],qfmsd=msd[4])
        all<-c(allmo,allsd)
      }else{
        all<-c(l4=lmo[1],etl4=lmo[2],rl4=lmo[3],ql4=lmo[4],
               fm=mo[1],etfm=mo[2],rfm=mo[3],qfm=mo[4])
      }
      return(all)
    }
    
    NRSssimple<-function (x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray"),sd=FALSE){
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
      
      mmm1<-mmm(x=sortedx,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
      rqscale1<-rqscale(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmscale,drm=drmscale,dlqm=dlqmscale,dqm=dqmscale,sd=sd)
      rqtm1<-rqtm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmtm,drm=drmtm,dlqm=dlqmtm,dqm=dqmtm,sd=sd,drsd=drmscale,dqsd=dqmscale)
      rqfm1<-rqfm(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,dlrm=dlrmfm,drm=drmfm,dlqm=dlqmfm,dqm=dqmfm,sd=sd,drsd=drmscale,dqsd=dqmscale)
      
      if(sd){
        first<-c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4])
        second<-c(l2=rqscale1[1],etl2=rqscale1[2],rl2=rqscale1[3],ql2=rqscale1[4],sd=sqrt(rqscale1[5]),etsd=sqrt(rqscale1[6]),
                  rsd=sqrt(rqscale1[7]),qsd=sqrt(rqscale1[8]))
        third<-c(l3=rqtm1[1]/rqscale1[1],etl3=rqtm1[2]/rqscale1[2],rl3=rqtm1[3]/rqscale1[3],ql3=rqtm1[4]/rqscale1[4],
                 skew=(rqtm1[5])/((rqscale1[5])^(3/2)),etskew=(rqtm1[6])/((rqscale1[6])^(3/2)),
                 rskew=(rqtm1[7])/((rqscale1[7])^(3/2)),qskew=(rqtm1[8])/((rqscale1[8])^(3/2)))
        fourth<-c(l4=rqfm1[1]/rqscale1[1],etl4=rqfm1[2]/rqscale1[2],rl4=rqfm1[3]/rqscale1[3],ql4=rqfm1[4]/rqscale1[4],kurt=(rqfm1[5])/((rqscale1[5])^(2)),etkurt=(rqfm1[6])/((rqscale1[6])^(2)),rkurt=(rqfm1[7])/((rqscale1[7])^(2)),qkurt=(rqfm1[8])/((rqscale1[8])^(2)))
        
        allmo<-list(first=first,second=second,third=third,fourth=fourth)
        
        
        firstsd<-c(sd=sqrt(rqscale1[5]),etsd=sqrt(rqscale1[6]),rsd=sqrt(rqscale1[7]),qsd=sqrt(rqscale1[8]))
        secondsd<-c(l2sd=rqscale1[9],etl2sd=rqscale1[10],rl2sd=rqscale1[11],ql2sd=rqscale1[12],
                    sdsd=second[5]*(1/2)*(rqscale1[13]/rqscale1[5]),
                    etsdsd=second[6]*(1/2)*(rqscale1[14]/rqscale1[6]),rsdsd=second[7]*(1/2)*(rqscale1[15]/rqscale1[7]),qsdsd=second[8]*(1/2)*(rqscale1[16]/rqscale1[8]))
        
        thirdsd<-c(l3sd=(rqtm1[1]/rqscale1[1])*((rqtm1[9]/rqtm1[1])^2+(rqscale1[9]/rqscale1[1])^2)^(1/2),
                   etl3sd=(rqtm1[2]/rqscale1[2])*((rqtm1[10]/rqtm1[2])^2+(rqscale1[10]/rqscale1[2])^2)^(1/2),
                   rl3sd=(rqtm1[3]/rqscale1[3])*((rqtm1[11]/rqtm1[3])^2+(rqscale1[11]/rqscale1[3])^2)^(1/2),
                   ql3sd=(rqtm1[4]/rqscale1[4])*((rqtm1[12]/rqtm1[4])^2+(rqscale1[12]/rqscale1[4])^2)^(1/2),
                   skewsd=third[5]*((rqtm1[13]/rqtm1[5])^2+(((((rqscale1[5])^(3/2))*(3/2)*(rqscale1[13]/rqscale1[5]))/((rqscale1[5])^(3/2))))^2)^(1/2),
                   etskewsd=third[6]*((rqtm1[14]/rqtm1[6])^2+(((((rqscale1[6])^(3/2))*(3/2)*(rqscale1[14]/rqscale1[6]))/((rqscale1[6])^(3/2))))^2)^(1/2),
                   rskewsd=third[7]*((rqtm1[15]/rqtm1[7])^2+(((((rqscale1[7])^(3/2))*(3/2)*(rqscale1[15]/rqscale1[7]))/((rqscale1[7])^(3/2))))^2)^(1/2),
                   qskewsd=third[8]*((rqtm1[16]/rqtm1[8])^2+(((((rqscale1[8])^(3/2))*(3/2)*(rqscale1[16]/rqscale1[8]))/((rqscale1[8])^(3/2))))^2)^(1/2))
        
        fourthsd<-c(l4sd=(rqfm1[1]/rqscale1[1])*((rqfm1[9]/rqfm1[1])^2+(rqscale1[9]/rqscale1[1])^2)^(1/2),
                    etl4sd=(rqfm1[2]/rqscale1[2])*((rqfm1[10]/rqfm1[2])^2+(rqscale1[10]/rqscale1[2])^2)^(1/2),
                    rl4sd=(rqfm1[3]/rqscale1[3])*((rqfm1[11]/rqfm1[3])^2+(rqscale1[11]/rqscale1[3])^2)^(1/2),
                    ql4sd=(rqfm1[4]/rqscale1[4])*((rqfm1[12]/rqfm1[4])^2+(rqscale1[12]/rqscale1[4])^2)^(1/2),
                    kurtsd=fourth[5]*((rqfm1[13]/rqfm1[5])^2+(((((rqscale1[5])^(2))*(2)*(rqscale1[13]/rqscale1[5]))/((rqscale1[5])^(2))))^2)^(1/2),
                    etkurtsd=fourth[6]*((rqfm1[14]/rqfm1[6])^2+(((((rqscale1[6])^(2))*(2)*(rqscale1[14]/rqscale1[6]))/((rqscale1[6])^(2))))^2)^(1/2),
                    rkurtsd=fourth[7]*((rqfm1[15]/rqfm1[7])^2+(((((rqscale1[7])^(2))*(2)*(rqscale1[15]/rqscale1[7]))/((rqscale1[7])^(2))))^2)^(1/2),
                    qkurtsd=fourth[8]*((rqfm1[16]/rqfm1[8])^2+(((((rqscale1[8])^(2))*(2)*(rqscale1[16]/rqscale1[8]))/((rqscale1[8])^(2))))^2)^(1/2))
        
        allsd<-list(firstsd=firstsd,secondsd=secondsd,thirdsd=thirdsd,fourthsd=fourthsd)
        
        all<-c(allmo,allsd)
      }else{
        first<-c(mean=mmm1[1],etm=mmm1[2],rm=mmm1[3],qm=mmm1[4])
        second<-c(l2=rqscale1[1],etl2=rqscale1[2],rl2=rqscale1[3],ql2=rqscale1[4],sd=sqrt(rqscale1[5]),etsd=sqrt(rqscale1[6]),
                  rsd=sqrt(rqscale1[7]),qsd=sqrt(rqscale1[8]))
        third<-c(l3=rqtm1[1]/rqscale1[1],etl3=rqtm1[2]/rqscale1[2],rl3=rqtm1[3]/rqscale1[3],ql3=rqtm1[4]/rqscale1[4],
                 skew=(rqtm1[5])/((rqscale1[5])^(3/2)),etskew=(rqtm1[6])/((rqscale1[6])^(3/2)),
                 rskew=(rqtm1[7])/((rqscale1[7])^(3/2)),qskew=(rqtm1[8])/((rqscale1[8])^(3/2)))
        fourth<-c(l4=rqfm1[1]/rqscale1[1],etl4=rqfm1[2]/rqscale1[2],rl4=rqfm1[3]/rqscale1[3],ql4=rqfm1[4]/rqscale1[4],kurt=(rqfm1[5])/((rqscale1[5])^(2)),etkurt=(rqfm1[6])/((rqscale1[6])^(2)),rkurt=(rqfm1[7])/((rqscale1[7])^(2)),qkurt=(rqfm1[8])/((rqscale1[8])^(2)))
        all<-list(first=first,second=second,third=third,fourth=fourth)
      }
      
      if((rqfm1[8])/((rqscale1[8])^(2))<4.7 & standist=="exponential"){
        print("The quantile kurtosis is lower than 4.7, it might be better to use the Rayleigh distribution as the standard distribution.")
      }
      return(all)
    }
    effectsizeNRSssimple<-function(x,y,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray")){
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
      lengthx<-length(sortedx)
      sortedy<-sort(y,decreasing = FALSE,method ="radix")
      lengthy<-length(sortedy)
      estimatex<-NRSssimple(x=sortedx,interval=interval,fast=fast,batch=batch,boot=boot,times =times,standist=standist,sd=TRUE)
      estimatey<-NRSssimple(x=sortedy,interval=interval,fast=fast,batch=batch,boot=boot,times =times,standist=standist,sd=TRUE)
      firsteffectsize<-c(mean=(estimatex$first[1]-estimatey$first[1])/(((((estimatex$first[1])^2)+((estimatey$first[1])^2))*0.5)^0.5),etm=(estimatex$first[2]-estimatey$first[2])/(((((estimatex$first[2])^2)+((estimatey$first[2])^2))*0.5)^0.5),rm=(estimatex$first[3]-estimatey$first[3])/(((((estimatex$first[3])^2)+((estimatey$first[3])^2))*0.5)^0.5),qm=(estimatex$first[4]-estimatey$first[4])/(((((estimatex$first[4])^2)+((estimatey$first[4])^2))*0.5)^0.5))
      
      secondeffectsize<-c(l2=(estimatex$second[1]-estimatey$second[1])/(((((estimatex$second[1])^2)+((estimatey$second[1])^2))*0.5)^0.5),etl2=(estimatex$second[2]-estimatey$second[2])/(((((estimatex$second[2])^2)+((estimatey$second[2])^2))*0.5)^0.5),rl2=(estimatex$second[3]-estimatey$second[3])/(((((estimatex$second[3])^2)+((estimatey$second[3])^2))*0.5)^0.5),ql2=(estimatex$second[4]-estimatey$second[4])/(((((estimatex$second[4])^2)+((estimatey$second[4])^2))*0.5)^0.5),
                          sd=(estimatex$second[5]-estimatey$second[5])/(((((estimatex$second[5])^2)+((estimatey$second[5])^2))*0.5)^0.5),etsd=(estimatex$second[6]-estimatey$second[6])/(((((estimatex$second[6])^2)+((estimatey$second[6])^2))*0.5)^0.5),rsd=(estimatex$second[7]-estimatey$second[7])/(((((estimatex$second[7])^2)+((estimatey$second[7])^2))*0.5)^0.5),qsd=(estimatex$second[8]-estimatey$second[8])/(((((estimatex$second[8])^2)+((estimatey$second[8])^2))*0.5)^0.5)
      )
      thirdeffectsize<-c(l3=(estimatex$third[1]-estimatey$third[1])/(((((estimatex$third[1])^2)+((estimatey$third[1])^2))*0.5)^0.5),etl3=(estimatex$third[2]-estimatey$third[2])/(((((estimatex$third[2])^2)+((estimatey$third[2])^2))*0.5)^0.5),rl3=(estimatex$third[3]-estimatey$third[3])/(((((estimatex$third[3])^2)+((estimatey$third[3])^2))*0.5)^0.5),ql3=(estimatex$third[4]-estimatey$third[4])/(((((estimatex$third[4])^2)+((estimatey$third[4])^2))*0.5)^0.5),
                         skew=(estimatex$third[5]-estimatey$third[5])/(((((estimatex$third[5])^2)+((estimatey$third[5])^2))*0.5)^0.5),etskew=(estimatex$third[6]-estimatey$third[6])/(((((estimatex$third[6])^2)+((estimatey$third[6])^2))*0.5)^0.5),rskew=(estimatex$third[7]-estimatey$third[7])/(((((estimatex$third[7])^2)+((estimatey$third[7])^2))*0.5)^0.5),qskew=(estimatex$third[8]-estimatey$third[8])/(((((estimatex$third[8])^2)+((estimatey$third[8])^2))*0.5)^0.5)
      )
      fourtheffectsize<-c(l4=(estimatex$fourth[1]-estimatey$fourth[1])/(((((estimatex$fourth[1])^2)+((estimatey$fourth[1])^2))*0.5)^0.5),etl4=(estimatex$fourth[2]-estimatey$fourth[2])/(((((estimatex$fourth[2])^2)+((estimatey$fourth[2])^2))*0.5)^0.5),rl4=(estimatex$fourth[3]-estimatey$fourth[3])/(((((estimatex$fourth[3])^2)+((estimatey$fourth[3])^2))*0.5)^0.5),ql4=(estimatex$fourth[4]-estimatey$fourth[4])/(((((estimatex$fourth[4])^2)+((estimatey$fourth[4])^2))*0.5)^0.5),
                          kurt=(estimatex$fourth[5]-estimatey$fourth[5])/(((((estimatex$fourth[5])^2)+((estimatey$fourth[5])^2))*0.5)^0.5),etkurt=(estimatex$fourth[6]-estimatey$fourth[6])/(((((estimatex$fourth[6])^2)+((estimatey$fourth[6])^2))*0.5)^0.5),rkurt=(estimatex$fourth[7]-estimatey$fourth[7])/(((((estimatex$fourth[7])^2)+((estimatey$fourth[7])^2))*0.5)^0.5),qkurt=(estimatex$fourth[8]-estimatey$fourth[8])/(((((estimatex$fourth[8])^2)+((estimatey$fourth[8])^2))*0.5)^0.5)
      )
      
      all<-list(firsteffectsize=firsteffectsize,secondeffectsize=secondeffectsize,thirdeffectsize=thirdeffectsize,fourtheffectsize=fourtheffectsize,estimatex=estimatex,estimatey=estimatey)
      if(estimatex$fourth[5]<4.7 & standist=="exponential"){
        print("The quantile kurtosis is lower than 4.7, it might be better to use the Rayleigh distribution as the standard distribution.")
      }
      
      return(all)
    }
    
    effectsizeNRSs1<-effectsizeNRSssimple(x=datax[i,],y=datay[i,],interval=interval,fast=fast,batch=batch,boot=boot,times =times,standist=standist)
    estimate1<-c(effectsizeNRSs1$firsteffectsize,effectsizeNRSs1$secondeffectsize,effectsizeNRSs1$thirdeffectsize,effectsizeNRSs1$fourtheffectsize)
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
  
  bootlist1a<-sort(bootlist1)
  bootlist2a<-sort(bootlist2)
  bootlist3a<-sort(bootlist3)
  bootlist4a<-sort(bootlist4)
  bootlist5a<-sort(bootlist5)
  bootlist6a<-sort(bootlist6)
  bootlist7a<-sort(bootlist7)
  bootlist8a<-sort(bootlist8)
  bootlist9a<-sort(bootlist9)
  bootlist10a<-sort(bootlist10)
  bootlist11a<-sort(bootlist11)
  bootlist12a<-sort(bootlist12)
  bootlist13a<-sort(bootlist13)
  bootlist14a<-sort(bootlist14)
  bootlist15a<-sort(bootlist15)
  bootlist16a<-sort(bootlist16)
  bootlist17a<-sort(bootlist17)
  bootlist18a<-sort(bootlist18)
  bootlist19a<-sort(bootlist19)
  bootlist20a<-sort(bootlist20)
  bootlist21a<-sort(bootlist21)
  bootlist22a<-sort(bootlist22)
  
  bootlist23a<-sort(bootlist23)
  bootlist24a<-sort(bootlist24)
  bootlist25a<-sort(bootlist25)
  bootlist26a<-sort(bootlist26)
  bootlist27a<-sort(bootlist27)
  bootlist28a<-sort(bootlist28)

  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  effectsizeNRSs2<-effectsizeNRSssimple(x=sortedx,y=sortedy,interval=interval,fast=fast,batch=batch,boot=boot,times =times,standist=standist)
  estimate<-c(effectsizeNRSs2$firsteffectsize,effectsizeNRSs2$secondeffectsize,effectsizeNRSs2$thirdeffectsize,effectsizeNRSs2$fourtheffectsize)
  
  bootlist1a<-bootlist1a-(estimate[1])
  bootlist2a<-bootlist2a-(estimate[2])
  bootlist3a<-bootlist3a-(estimate[3])
  bootlist4a<-bootlist4a-(estimate[4])
  bootlist5a<-bootlist5a-(estimate[5])
  bootlist6a<-bootlist6a-(estimate[6])
  bootlist7a<-bootlist7a-(estimate[7])
  bootlist8a<-bootlist8a-(estimate[8])
  bootlist9a<-bootlist9a-(estimate[9])
  bootlist10a<-bootlist10a-(estimate[10])
  bootlist11a<-bootlist11a-(estimate[11])
  bootlist12a<-bootlist12a-(estimate[12])
  bootlist13a<-bootlist13a-(estimate[13])
  bootlist14a<-bootlist14a-(estimate[14])
  bootlist15a<-bootlist15a-(estimate[15])
  bootlist16a<-bootlist16a-(estimate[16])
  bootlist17a<-bootlist17a-(estimate[17])
  bootlist18a<-bootlist18a-(estimate[18])
  bootlist19a<-bootlist19a-(estimate[19])
  bootlist20a<-bootlist20a-(estimate[20])
  bootlist21a<-bootlist21a-(estimate[21])
  bootlist22a<-bootlist22a-(estimate[22])
  bootlist23a<-bootlist23a-(estimate[23])
  bootlist24a<-bootlist24a-(estimate[24])
  bootlist25a<-bootlist25a-(estimate[25])
  bootlist26a<-bootlist26a-(estimate[26])
  bootlist27a<-bootlist27a-(estimate[27])
  bootlist28a<-bootlist28a-(estimate[28])
  
  cidiff<-function(bootlist,low,up,estimate,n){
    (estimate[n])+c(bootlist[low],bootlist[up])
  }
  
  ci_diff<-c(c(mean=cidiff(bootlist1a,low=low,up=up,estimate=estimate,n=1),etm=cidiff(bootlist2a,low=low,up=up,estimate=estimate,n=2),rm=cidiff(bootlist3a,low=low,up=up,estimate=estimate,n=3),qm=cidiff(bootlist4a,low=low,up=up,estimate=estimate,n=4)),
             c(l2=cidiff(bootlist5a,low=low,up=up,estimate=estimate,n=5),etl2=cidiff(bootlist6a,low=low,up=up,estimate=estimate,n=6),rl2=cidiff(bootlist7a,low=low,up=up,estimate=estimate,n=7),ql2=cidiff(bootlist8a,low=low,up=up,estimate=estimate,n=8),sd=cidiff(bootlist9a,low=low,up=up,estimate=estimate,n=9),etsd=cidiff(bootlist10a,low=low,up=up,estimate=estimate,n=10),
               rsd=cidiff(bootlist11a,low=low,up=up,estimate=estimate,n=11),qsd=cidiff(bootlist12a,low=low,up=up,estimate=estimate,n=12)),
             c(l3=cidiff(bootlist13a,low=low,up=up,estimate=estimate,n=13),etl3=cidiff(bootlist14a,low=low,up=up,estimate=estimate,n=14),rl3=cidiff(bootlist15a,low=low,up=up,estimate=estimate,n=15),ql3=cidiff(bootlist16a,low=low,up=up,estimate=estimate,n=16),
               skew=cidiff(bootlist17a,low=low,up=up,estimate=estimate,n=17),etskew=cidiff(bootlist18a,low=low,up=up,estimate=estimate,n=18),
               rskew=cidiff(bootlist19a,low=low,up=up,estimate=estimate,n=19),qskew=cidiff(bootlist20a,low=low,up=up,estimate=estimate,n=20)),
             c(l4=cidiff(bootlist21a,low=low,up=up,estimate=estimate,n=21),etl4=cidiff(bootlist22a,low=low,up=up,estimate=estimate,n=22),rl4=cidiff(bootlist23a,low=low,up=up,estimate=estimate,n=23),ql4=cidiff(bootlist24a,low=low,up=up,estimate=estimate,n=24),kurt=cidiff(bootlist25a,low=low,up=up,estimate=estimate,n=25),etkurt=cidiff(bootlist26a,low=low,up=up,estimate=estimate,n=26),rkurt=cidiff(bootlist27a,low=low,up=up,estimate=estimate,n=27),qkurt=cidiff(bootlist28a,low=low,up=up,estimate=estimate,n=28)))
  
  
  se<-c(c(mean=sd(bootlist1a),etm=sd(bootlist2a),rm=sd(bootlist3a),qm=sd(bootlist4a)),
             c(l2=sd(bootlist5a),etl2=sd(bootlist6a),rl2=sd(bootlist7a),ql2=sd(bootlist8a),sd=sd(bootlist9a),etsd=sd(bootlist10a),
               rsd=sd(bootlist11a),qsd=sd(bootlist12a)),
             c(l3=sd(bootlist13a),etl3=sd(bootlist14a),rl3=sd(bootlist15a),ql3=sd(bootlist16a),
               skew=sd(bootlist17a),etskew=sd(bootlist18a),
               rskew=sd(bootlist19a),qskew=sd(bootlist20a)),
             c(l4=sd(bootlist21a),etl4=sd(bootlist22a),rl4=sd(bootlist23a),ql4=sd(bootlist24a),kurt=sd(bootlist25a),etkurt=sd(bootlist26a),rkurt=sd(bootlist27a),qkurt=sd(bootlist28a)))
  all<-list(ci_diff=ci_diff,se=se,estimate=estimate)
  
  return(all)
}

effectsizeNRSs<-function(x,y,ci=TRUE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist=c("exponential","Rayleigh","exp","Ray"),alpha=0.05,nboot=100){
  if (ci){
    return(esbootparallel(x=x,y=y,interval=interval,fast=fast,batch=batch,boot=boot,times =times,standist=standist,alpha=alpha,nboot=nboot))
  } 
  else{return (effectsizeNRSssimple(x=x,y=y,interval=interval,fast=fast,batch=batch,boot=boot,times =times,standist=standist))
  }}
winsor<-function (x, fraction=1/9)
{
  lim <- quantile(x, probs=c(fraction, 1-fraction))
  x[ x < lim[1] ] <- lim[1]
  x[ x > lim[2] ] <- lim[2]
  x
}
twmean1<-function(x,fraction=1/9,type=1){
  c(mean=mean(x),trimmedmean=mean(x,fraction),winsorizedmean=mean(winsor(x,fraction)))[type]
}
twmean<-function(x,fraction=1/9){
  c(mean=mean(x),trimmedmean=mean(x,fraction),winsorizedmean=mean(winsor(x,fraction)))
}
twcov<-function (x,y=NULL){
  matrix1 <- cbind(x, y)
  matrix1 <- na.omit(matrix1)
  x <- matrix1[, 1]
  y <- matrix1[, 2]
  productxy<-x*y
  productx2<-x*x
  producty2<-y*y
  meanx <- twmean(x,fraction=1/9)
  meany <- twmean(y,fraction=1/9)
  meanxy <- twmean(productxy,fraction=1/9)
  meanx2 <- twmean(productx2,fraction=1/9)
  meany2 <- twmean(producty2,fraction=1/9)
  mean1<-(meanxy[1]-meanx[1]*meany[1])
  tm1<-(meanxy[2]-meanx[2]*meany[2])
  wm1<-(meanxy[3]-meanx[3]*meany[3])
  twcov <- c(mean1=mean1,tm1=tm1,wm1=wm1)
  return(twcov)
}
twreg<-function(x,y=NULL,iter = 20){
  x <- as.matrix(x)
  xy <- cbind(x, y)
  xy <- na.omit(xy)
  ncolx<-ncol(x)
  ncolx1<-ncol(x)+1
  x <- xy[, 1:ncolx]
  y <- xy[, ncolx1]
  x = as.matrix(x)
  matrix1 <- matrix(0, ncol(x), 1)
  matrix2 <- matrix(0, ncol(x), ncol(x))
  coef<-c()
  for (type in (1:3)){
    mval1 <- apply(x, 2, twmean1,fraction=1/9,type=type)
    for (i in 1:ncol(x)) {
      matrix1[i, 1] <- twcov(x=x[, i],y=y)[type]
      for (j in 1:ncol(x)) matrix2[i, j] <-twcov(x=x[, i],y=x[, j])[type]
    }
    slope <- solve(matrix2, matrix1)
    inter <- twmean1(x=y,fraction=1/9,type=type) - sum(slope %*% mval1)
    for (it in 1:iter) {
      res <- y - x %*% slope - inter
      for (i in 1:ncol(x)) matrix1[i, 1] <- twcov(x=x[, i],y=res)[type]
      sadd <- solve(matrix2, matrix1)
      inadd <- twmean1(x=res,fraction=1/9,type=type) - sum(sadd %*% mval1)
      if (max(abs(sadd), abs(inadd)) < 1e-04) 
        break
      slope <- slope + sadd
      inter <- inter + inadd
    }
    if (max(abs(sadd), abs(inadd)) >= 1e-04) {
      namestype<-c("mean", "tm", "wm")  
      print(paste(namestype[type],"failed to converge within", iter, "iterations"))}
    all1<-(c(Intercept=inter, slope=slope,ResidualSE=sd(res)))
    coef<-rbind(coef,all1)
  }
  rownames(coef) <- c("mean", "tm", "wm")  
  coef
}



rqcov<-function (x,y=NULL,interval=9,fast=TRUE,batch="auto",standist=c("exponential","Rayleigh","exp","Ray")){
  matrix1 <- cbind(x, y)
  matrix1 <- na.omit(matrix1)
  x <- matrix1[, 1]
  y <- matrix1[, 2]
  productxy<-x*y
  productx2<-x*x
  producty2<-y*y
  if(standist=="exponential"|| standist=="exp"){
    drm=0.3665
    dqm=0.82224

  }else if (standist=="Rayleigh"|| standist=="Ray"){
    drm=0.4025526
    dqm=0.4452798
  }
  meanx <- rqmean(x,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  meany <- rqmean(y,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  meanxy <- rqmean(productxy,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  meanx2 <- rqmean(productx2,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  meany2 <- rqmean(producty2,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm)
  mean1<-(meanxy[1]-meanx[1]*meany[1])
  etm1<-(meanxy[2]-meanx[2]*meany[2])
  rm1<-(meanxy[3]-meanx[3]*meany[3])
  qm1<-(meanxy[4]-meanx[4]*meany[4])
  rqcov <- c(mean1=mean1,etm1=etm1,rm1=rm1,qm1=qm1)
  return(rqcov)
}
rqreg<-function(x,y=NULL,iter = 20,interval=9,fast=TRUE,batch="auto",standist=c("exponential","Rayleigh","exp","Ray")){
  if(standist=="exponential"|| standist=="exp"){
    drm=0.3665
    dqm=0.82224
    
  }else if (standist=="Rayleigh"|| standist=="Ray"){
    drm=0.4025526
    dqm=0.4452798
  }
  x <- as.matrix(x)
  xy <- cbind(x, y)
  xy <- na.omit(xy)
  ncolx<-ncol(x)
  ncolx1<-ncol(x)+1
  x <- xy[, 1:ncolx]
  y <- xy[, ncolx1]
  x = as.matrix(x)
  matrix1 <- matrix(0, ncol(x), 1)
  matrix2 <- matrix(0, ncol(x), ncol(x))
  coef<-c()
  for (type in (1:4)){
    mval1 <- apply(x, 2, mmme,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm,type=type)
    for (i in 1:ncol(x)) {
      matrix1[i, 1] <- rqcov(x=x[, i],y=y,interval=interval,fast=fast,batch=batch,standist=standist)[type]
      for (j in 1:ncol(x)) matrix2[i, j] <-rqcov(x=x[, i],y=x[, j],interval=interval,fast=fast,batch=batch,standist=standist)[type]
    }
    slope <- solve(matrix2, matrix1)
    inter <- mmme(x=y,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm,type=type) - sum(slope %*% mval1)
    for (it in 1:iter) {
      res <- y - x %*% slope - inter
      for (i in 1:ncol(x)) matrix1[i, 1] <- rqcov(x=x[, i],y=res,interval=interval,fast=fast,batch=batch,standist=standist)[type]
      sadd <- solve(matrix2, matrix1)
      inadd <- mmme(x=res,interval=interval,fast=fast,batch=batch,drm=drm,dqm=dqm,type=type) - sum(sadd %*% mval1)
      if (max(abs(sadd), abs(inadd)) < 1e-04) 
        break
      slope <- slope + sadd
      inter <- inter + inadd
    }
    if (max(abs(sadd), abs(inadd)) >= 1e-04) {
      namestype<-c("mean", "etm", "rm", "qm")  
      print(paste(namestype[type],"failed to converge within", iter, "iterations"))}
    all1<-(c(Intercept=inter, slope=slope,ResidualSE=sd(res)))
    coef<-rbind(coef,all1)
  }
  rownames(coef) <- c("mean", "etm", "rm", "qm")  
  coef
}
library(lmtest)
#robust regression test
#Gaussian outliers
n<-5400
y<-as.numeric(n)
x<-as.numeric(n)
error<-as.numeric(n)
for (i in 1:n){
  x1 <- rnorm(1,0,1)
  x2 <- runif(1,200,201)
  u <- runif(1)
  k <- as.integer(u > 0.99) 
  error[i] <- (1-k)* x1 +  k* x2 
  x[i]<-runif(1,0,10)
  y[i]<-10+2*x[i]+error[i]
}
hist(error)
m1=lm(y~x)
summary(m1)
rqreg(x=x, y=y,iter = 200,interval=9,fast=TRUE,batch="auto",standist="exp")
twreg(x=x, y=y,iter = 200)


#exponential outliers
n<-5400
y<-as.numeric(n)
x<-as.numeric(n)
error<-as.numeric(n)
for (i in 1:n){
  x1 <- rexp(1,1)-1
  x2 <- runif(1,200,201)
  u <- runif(1)
  k <- as.integer(u > 0.99) 
  error[i] <- (1-k)* x1 +  k* x2 
  x[i]<-runif(1,0,10)
  y[i]<-10+2*x[i]+error[i]
}
hist(error)
m1=lm(y~x)
summary(m1)
rqreg(x=x, y=y,iter = 200,interval=9,fast=TRUE,batch="auto",standist="exp")
twreg(x=x, y=y,iter = 200)

#Rayleigh outliers
n<-5400
y<-as.numeric(n)
x<-as.numeric(n)
error<-as.numeric(n)
for (i in 1:n){
  x1 <- rRayleigh(1,1)-sqrt(pi/2)
  x2 <- runif(1,200,201)
  u <- runif(1)
  k <- as.integer(u > 0.99) 
  error[i] <- (1-k)* x1 +  k* x2 
  x[i]<-runif(1,0,10)
  y[i]<-10+2*x[i]+error[i]
}
hist(error)
m1=lm(y~x)
summary(m1)
rqreg(x=x, y=y,iter = 200,interval=9,fast=TRUE,batch="auto",standist="exp")
twreg(x=x, y=y,iter = 200)
#the performance of rm and qm is not as good as etm


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

#the bootstrap method for confidential interval, standard errors and hypothesis testing, are from three textbooks, including 
#Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.
#Bickel, P. J., & Doksum, K. A. (2015). Mathematical statistics: basic ideas and selected topics, volumes I-II package. Chapman and Hall/CRC.
#Rice, J. A. (2006). Mathematical statistics and data analysis. Cengage Learning.

#the asymptotic validity of bootstrap of U-statistics, although haven't been proven yet, is verified by Monta Carlo study and highly suggested by 
#Bickel, P. J., & Freedman, D. A. (1981). Some asymptotic theory for the bootstrap. The annals of statistics, 9(6), 1196-1217.
#Bickel, P. J., & Freedman, D. A. (1984). Asymptotic normality and the bootstrap in stratified sampling. The annals of statistics, 470-482.

#this standard deviation of the distribution of U-statistic is calculated based on the law of prorogation of uncertainty.
NRSs(x=xexp,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exp",cise = FALSE,parallel=TRUE,alpha = 0.05,nboot = 100, sd=TRUE)
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

#To make comparisons easier, sample standardized moments and scaled L-moments are provided (this moments are estimated same based on U-statistics, not the formula approach
#because if direct compare sample moments with NRSs, the variance from bootstrap also need to take into consideration.

#ETM-based moments and L-moments are also provided, although the biases are large compared to NRSs, ET-moments have near-optimum standard errors and so can be expected will have a place in hypothesis testing.

#The standard deviations of the distributions of U-statistic can be used to estimate the consistency percentage.

#another very interesting and probably novel approach is effect size based on the standard deviation of the U-statistics (I said novel is because I haven't found effect size of higher moments, e.g., standard deviation, skewness, kurtosis)
xexp<-c(rexp(5400,2))
yexp<-rexp(5400,1)

effectsizeNRSs(x=xexp,y=yexp,ci=FALSE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exp",alpha=0.05,nboot=100)
#the standard deviation is estimated by the corresponding estimators, e.g., for rskew, the standard deviation of the distribution of U-statistics is estimated by rsd.

#the confidence interval of effect size can be estimated by bootstrap. With the help of parallel computing, the running time is around 1 mins
#because for every estimators, the standard deviation of the corresponding U statistics needs to be estimated by etsd, rsd, qsd.

#to reduce the test time, the boot times of U-statistics are 5400, instead of 54000.
effectsizeNRSs(x=xexp,y=yexp,ci=TRUE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)
#the se is the standard error of effect size.

#The standard error and confidential interval of the robust or quantile mean can be accurately estimated by bootstrapping.
xexp<-c(rexp(5400,1))
rqmean(x=xexp,interval=9,fast=TRUE,batch="auto",drm=0.3665,dqm=0.82224,cise = TRUE,alpha = 0.05,nboot = 1000)
#A similar approach can be applied to all NRSs, but just 100 nboot takes ~10 mins.

#If you don't want to wait for a long time to compare, don't run the following code.
#NRSs(x,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,parallel=FALSE,alpha = 0.05,nboot = 100)

#A solution is parallel computing (takes 1 min with 12 cores, but is unavailable on some types of computers).

#The standard errors of robust skewness and kurtosis are lower than those of sample skewness and kurtosis.

#also, the one-sample hypothesis testing can be done with the equal-tail bootstrap P-value.

#the null values are the corresponding population parameters.

#the P value is two-side.

NRSs(x=xexp,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,parallel=TRUE,alpha = 0.05,nboot=100,null_mean=1,null_sd=1,null_skew=2,null_kurt=9,null_l2=0.5,null_l3=1/3,null_l4=1/6)


#two-goup comparison can also be done with a similar approach.
xexp<-rexp(5400,1.1)
yexp<-rexp(5400,1)

#to reduce the test time, the boot times of U-statistics are 5400, instead of 54000.
#empirical bootstrap hypothesis test
htest(x=xexp,y=yexp,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)
#percentile bootstrap hypothesis test (very controversial, but the results are similar.)
htest(x=xexp,y=yexp,boottype="percentile",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

xexp<-c(rexp(5350,1),rnorm(50,30))
yexp<-rexp(5400,1)

#test of outliers.

htest(x=xexp,y=yexp,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

htest(x=xexp,y=yexp,boottype="percentile",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

xexp<-c(rexp(5400,1))
yexp<-rexp(5400,1)

#test of null hypothesis

htest(x=xexp,y=yexp,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)
htest(x=xexp,y=yexp,boottype="percentile",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)


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
NRSs(x=xgamma,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exponential",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=targetgam[1],null_sd=sqrt((a/100)),null_skew=2/sqrt(a/100),null_kurt=((6/(a/100))+3),null_l2=targetgam[2],null_l3=targetgam[3],null_l4=targetgam[4])

NRSs(x=xgamma,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=targetgam[1],null_sd=sqrt((a/100)),null_skew=2/sqrt(a/100),null_kurt=((6/(a/100))+3),null_l2=targetgam[2],null_l3=targetgam[3],null_l4=targetgam[4])

xgamma<-c(rgamma(5400, shape=a/100, rate = 2))
ygamma<-c(rgamma(5400, shape=a/100, rate = 1))

#test of effect size
effectsizeNRSs(x=xgamma,y=ygamma,ci=FALSE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

xgamma<-c(rgamma(5400, shape=a/100, rate = 1.1))
ygamma<-c(rgamma(5400, shape=a/100, rate = 1))

#empirical bootstrap hypothesis test
htest(x=xgamma,y=ygamma,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

xgamma<-c(rgamma(5400, shape=a/100, rate = 1))
ygamma<-c(rgamma(5400, shape=a/100, rate = 1))

#test of null hypothesis
htest(x=xgamma,y=ygamma,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

a=150
xgamma<-c(rgamma(5400, shape=a/100, rate = 1))
targetgam<-lmrgam(para = c(a/100, 1), nmom = 4)
NRSs(x=xgamma,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xgamma,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
#even the kurtosis is not very high, 7, the standard errors are still lower than sample moments. 
NRSs(x=xgamma,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=targetgam[1],null_sd=sqrt((a/100)),null_skew=2/sqrt(a/100),null_kurt=((6/(a/100))+3),null_l2=targetgam[2],null_l3=targetgam[3],null_l4=targetgam[4])

xgamma<-c(rgamma(5400, shape=a/100, rate = 2))
ygamma<-c(rgamma(5400, shape=a/100, rate = 1))

#test of effect size
effectsizeNRSs(x=xgamma,y=ygamma,ci=FALSE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

xgamma<-c(rgamma(5400, shape=a/100, rate = 1.1))
ygamma<-c(rgamma(5400, shape=a/100, rate = 1))

#empirical bootstrap hypothesis test
htest(x=xgamma,y=ygamma,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

xgamma<-c(rgamma(5400, shape=a/100, rate = 1))
ygamma<-c(rgamma(5400, shape=a/100, rate = 1))

#test of null hypothesis
htest(x=xgamma,y=ygamma,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)


xRayleigh<-rRayleigh(n=5400, scale = 1) 

NRSs(x=xRayleigh,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xRayleigh,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
NRSs(x=xRayleigh,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh",cise = TRUE,parallel=TRUE,alpha = 0.05,nboot=100,null_mean=sqrt(pi/2),null_sd=sqrt(2-(pi/2)),null_skew=2*sqrt((pi))*(pi-3)/((4-pi)^(3/2)),null_kurt=(3-(6*(pi)^2-24*(pi)+16)/((4-pi)^(2))),null_l2=0.5*(sqrt(2)-1)*sqrt(pi),null_l3=((1/6)*(2*sqrt(6)+3*sqrt(2)-9)*sqrt(pi))/(0.5*(sqrt(2)-1)*sqrt(pi)),null_l4=((sqrt((77/6)-5*sqrt(6))-3/4)*sqrt(2*pi))/(0.5*(sqrt(2)-1)*sqrt(pi)))


xRayleigh<-rRayleigh(n=5400, scale = 2) 
yRayleigh<-rRayleigh(n=5400, scale = 1) 

#test of effect size
effectsizeNRSs(x=xRayleigh,y=yRayleigh,ci=FALSE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

xRayleigh<-rRayleigh(n=5400, scale = 1.1) 
yRayleigh<-rRayleigh(n=5400, scale = 1) 

#empirical bootstrap hypothesis test
htest(x=xRayleigh,y=yRayleigh,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

xRayleigh<-rRayleigh(n=5400, scale = 1) 
yRayleigh<-rRayleigh(n=5400, scale = 1) 

#test of null hypothesis
htest(x=xRayleigh,y=yRayleigh,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

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



xpois<-rpois(5400,19)
ypois<-rpois(5400,8)

#test of effect size
effectsizeNRSs(x=xpois,y=ypois,ci=FALSE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

xpois<-rpois(5400,9)
ypois<-rpois(5400,8)

#empirical bootstrap hypothesis test
htest(x=xpois,y=ypois,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

xpois<-rpois(5400,8)
ypois<-rpois(5400,8)

#test of null hypothesis
htest(x=xpois,y=ypois,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)


xnorm<-c(rnorm(5400))

NRSs(x=xnorm,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xnorm,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")

NRSs(x=xnorm,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=0,null_sd=1,null_skew=0,null_kurt=3,null_l2=1/sqrt(pi),null_l3=0,null_l4=(30*(1/(pi))*(atan(sqrt(2)))-9))

xnorm<-c(rnorm(5400,2))
ynorm<-c(rnorm(5400,1))

#test of effect size
effectsizeNRSs(x=xnorm,y=ynorm,ci=FALSE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

xnorm<-c(rnorm(5400,1.1))
ynorm<-c(rnorm(5400,1))

#empirical bootstrap hypothesis test
htest(x=xnorm,y=ynorm,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

xnorm<-c(rnorm(5400,1))
ynorm<-c(rnorm(5400,1))

#test of null hypothesis
htest(x=xnorm,y=ynorm,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)


xlogis<-c(rlogis(5400, location = 0, scale = 1))

NRSs(x=xlogis,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xlogis,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
NRSs(x=xlogis,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=0,null_sd=sqrt(((pi^2)/3)),null_skew=0,null_kurt=(((6/5)+3)*((sqrt((pi^2)/3))^4))/(sqrt(((pi^2)/3))^(4)),null_l2=1,null_l3=0,null_l4=1/6)

xlogis<-c(rlogis(5400,2))
ylogis<-c(rlogis(5400,1))

#test of effect size
effectsizeNRSs(x=xlogis,y=ylogis,ci=FALSE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

xlogis<-c(rlogis(5400,1.1))
ylogis<-c(rlogis(5400,1))

#empirical bootstrap hypothesis test
htest(x=xlogis,y=ylogis,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

xlogis<-c(rlogis(5400,1))
ylogis<-c(rlogis(5400,1))

#test of null hypothesis
htest(x=xlogis,y=ylogis,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)


xlaplace<-c(rLaplace(n=5400, location = 0, scale = 1))
NRSs(x=xlaplace,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xlaplace,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
NRSs(x=xlaplace,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=0,null_sd=sqrt(2),null_skew=0,null_kurt=(6*(sqrt(2)^4))/(4),null_l2=3/4,null_l3=0,null_l4=1/(3*sqrt(2)))


xlaplace<-c(rLaplace(n=5400,location = 2,scale = 1))
ylaplace<-c(rLaplace(n=5400,location = 1,scale = 1))

#test of effect size
effectsizeNRSs(x=xlaplace,y=ylaplace,ci=FALSE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

xlaplace<-c(rLaplace(n=5400,location = 1.1,scale = 1))
ylaplace<-c(rLaplace(n=5400,location = 1,scale = 1))

#empirical bootstrap hypothesis test
htest(x=xlaplace,y=ylaplace,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)

xlaplace<-c(rLaplace(n=5400,location = 1,scale = 1))
ylaplace<-c(rLaplace(n=5400,location = 1,scale = 1))

#test of null hypothesis
htest(x=xlaplace,y=ylaplace,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="Rayleigh",alpha=0.05,nboot=100)


#NRSs have excellent performance even for heavy tailed distributions.

#two performance criteria, consistency (or sensitive) and standard error
library(lmom)
a=500
xpareto<-c(rPareto(5400, scale  = 1, shape=2+a/100))
targetlpareto<-lmrgpa(para = c(1,1/(2+a/100),- 1/(2+a/100)), nmom = 4)

NRSs(x=xpareto,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xpareto,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
#the standard errors are lower, especially for robust moments and L-moments
NRSs(x=xpareto,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=targetlpareto[1],null_sd=(sqrt(((2+a/100))*(1)/((-2+(2+a/100))*((-1+(2+a/100))^2)))),null_skew=((((2+a/100)+1)*(2)*(sqrt(a/100)))/((-3+(2+a/100))*(((2+a/100))^(1/2)))),null_kurt=((3+(6*((2+a/100)^3+(2+a/100)^2-6*(2+a/100)-2)/(((2+a/100))*((-3+(2+a/100)))*((-4+(2+a/100))))))),null_l2=targetlpareto[2],null_l3=targetlpareto[3],null_l4=targetlpareto[4])

xpareto<-c(rPareto(5400, scale  = 2, shape=2+a/100))
ypareto<-c(rPareto(5400, scale  = 1, shape=2+a/100))

#test of effect size
effectsizeNRSs(x=xpareto,y=ypareto,ci=FALSE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

xpareto<-c(rPareto(5400, scale  = 1.1, shape=2+a/100))
ypareto<-c(rPareto(5400, scale  = 1, shape=2+a/100))

#empirical bootstrap hypothesis test
htest(x=xpareto,y=ypareto,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

xpareto<-c(rPareto(5400, scale  = 1, shape=2+a/100))
ypareto<-c(rPareto(5400, scale  = 1, shape=2+a/100))

#test of null hypothesis
htest(x=xpareto,y=ypareto,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)


a=100
xlnorm<-c(rlnorm(5400,meanlog=0,sdlog=a/100))
targetlnorm<-lmrln3(para = c(0,0, a/100), nmom = 4)

NRSs(x=xlnorm,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xlnorm,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
#the standard errors are lower, especially for robust moments and L-moments
NRSs(x=xlnorm,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=targetlnorm[1],null_sd=sqrt((exp((a/100)^2)*(-1+exp((a/100)^2)))),null_skew=sqrt(exp((a/100)^2)-1)*((2+exp((a/100)^2))),null_kurt=(((-3+exp(4*((a/100)^2))+2*exp(3*((a/100)^2))+3*exp(2*((a/100)^2))))),null_l2=targetlnorm[2],null_l3=targetlnorm[3],null_l4=targetlnorm[4])

xlnorm<-c(rlnorm(5400,meanlog=2,sdlog=a/100))
ylnorm<-c(rlnorm(5400,meanlog=1,sdlog=a/100))

#test of effect size
effectsizeNRSs(x=xlnorm,y=ylnorm,ci=FALSE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

xlnorm<-c(rlnorm(5400,meanlog=1.1,sdlog=a/100))
ylnorm<-c(rlnorm(5400,meanlog=1,sdlog=a/100))

#empirical bootstrap hypothesis test
htest(x=xlnorm,y=ylnorm,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

xlnorm<-c(rlnorm(5400,meanlog=1,sdlog=a/100))
ylnorm<-c(rlnorm(5400,meanlog=1,sdlog=a/100))

#test of null hypothesis
htest(x=xlnorm,y=ylnorm,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)



a=150
xweibull<-c(rweibull(5400, shape=a/100, scale = 1))
library(lmom)
targetwei<-lmrwei(para = c(0, 1, a/100), nmom = 4)
NRSs(x=xweibull,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential")
NRSs(x=xweibull,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="Rayleigh")
NRSs(x=xweibull,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =54000,standist="exponential",cise = TRUE,alpha = 0.05,nboot = 100,null_mean=gamma(1+1/(a/100)),null_sd=(sqrt(gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2)),null_skew=(gamma(1+3/(a/100))-3*(gamma(1+1/(a/100)))*((gamma(1+2/(a/100))))+2*((gamma(1+1/(a/100)))^3))/((sqrt(gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2))^(3)),null_kurt=((gamma(1+4/(a/100))-4*(gamma(1+3/(a/100)))*((gamma(1+1/(a/100))))+6*(gamma(1+2/(a/100)))*((gamma(1+1/(a/100)))^2)-3*((gamma(1+1/(a/100)))^4))/(((gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2))^(2))),null_l2=targetwei[2],null_l3=targetwei[3],null_l4=targetwei[4])


xweibull<-c(rweibull(5400, shape=a/100, scale = 2))
yweibull<-c(rweibull(5400, shape=a/100, scale = 1))

#test of effect size
effectsizeNRSs(x=xweibull,y=yweibull,ci=FALSE,interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

xweibull<-c(rweibull(5400, shape=a/100, scale = 1.1))
yweibull<-c(rweibull(5400, shape=a/100, scale = 1))

#empirical bootstrap hypothesis test
htest(x=xweibull,y=yweibull,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)

xweibull<-c(rweibull(5400, shape=a/100, scale = 1))
yweibull<-c(rweibull(5400, shape=a/100, scale = 1))

#test of null hypothesis
htest(x=xweibull,y=yweibull,boottype="empirial",interval=9,fast=TRUE,batch="auto",boot=TRUE,times =5400,standist="exp",alpha=0.05,nboot=100)


#for more tests, use the codes in consistency.R
