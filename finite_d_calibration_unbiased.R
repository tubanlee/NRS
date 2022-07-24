batchnumber=1 


#it takes more than one day to run this file with default setting. Reduce the batchsizebase can reduce the time for easy testing
moments<-function (x){
  n<-length(x)
  m1<-mean(x)
  sd1<-sd(x)
  tm1<-(sum((x - m1)^3)/n)*(n^2/((n-1)*(n-2)))
  fm1<-(sum((x - m1)^4)/n)
  listall<-c(mean=m1,variance=(sd1)^2,tm=tm1,fm=fm1)
  (listall)
}
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

finddmmm<-function(expectboot,expecttrue,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  expectboot=mean(x)
  quatileexpectboot<-(min(which(sortedx>(expectboot)))-1)/length(x)
  etm1<-etm(x,interval=interval,fast=fast,batch=batch)
  mx1<-(min(which(sortedx>(etm1[2])))-1)/length(x)
  mx2<-1/2
  if (mx1>0.5){
    dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
  }else{
    quatileexpectboot<-1-quatileexpectboot
    mx1<-1-mx1
    dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
  }
  drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
  listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
  return(listd)
}

finddmmme<-function(expectboot,expecttrue,x,interval=9,fast=TRUE,batch=1000,sorted=FALSE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-sort(x,decreasing = FALSE,method ="radix")
  }
  lengthx<-length(x)
  quatileexpectboot<-(min(which(sortedx>(expectboot)))-1)/lengthx
  quatileexpecttrue<-(min(which(sortedx>(expecttrue)))-1)/lengthx
  etm1<-etm(x,interval=interval,fast=fast,batch=batch)
  mx1<-(min(which(sortedx>(etm1[2])))-1)/lengthx
  mx2<-1/2
  if (mx1>0.5){
    dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
  }else{
    quatileexpectboot<-1-quatileexpectboot
    quatileexpecttrue<-1-quatileexpecttrue
    mx1<-1-mx1
    dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
  }
  drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
  listd<-c(drm1,dqm1)
  return(listd)
}
finddtm<-function (x,expectbootl,expectboot,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE){
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
  dlmo<-finddmmm(expectboot=expectbootl,expecttrue=expectbootl,x=dp2lm,interval=9,fast=fast,batch=batch,sorted=FALSE)
  dmo<-finddmmm(expectboot=expectboot,expecttrue=expectboot,x=dp2m,interval=9,fast=fast,batch=batch,sorted=FALSE)
  all<-c(expectbootl=dlmo[1],etml=dlmo[2],ctml=dlmo[3],mx1l=dlmo[4],quatileexpectbootl=dlmo[5],
         expectboot=dmo[1],etm=dmo[2],ctm=dmo[3],mx1=dmo[4],quatileexpectboot=dmo[5]
  )
  return(all)
}

finddfm<-function (x,expectbootl,expectboot,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE){
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
  dlmo<-finddmmm(expectboot=expectbootl,expecttrue=expectbootl,x=dp2lm,interval=9,fast=fast,batch=batch,sorted=FALSE)
  dmo<-finddmmm(expectboot=expectboot,expecttrue=expectboot,x=dp2m,interval=9,fast=fast,batch=batch,sorted=FALSE)
  all<-c(expectbootl=dlmo[1],etml=dlmo[2],ctml=dlmo[3],mx1l=dlmo[4],quatileexpectbootl=dlmo[5],
         expectboot=dmo[1],etm=dmo[2],ctm=dmo[3],mx1=dmo[4],quatileexpectboot=dmo[5]
  )
  return(all)
}

finddscale<-function (x,expectbootl,expectboot,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE){
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
  dlmo<-finddmmm(expectboot=expectbootl,expecttrue=expectbootl,x=dp2lm,interval=9,fast=TRUE,batch=10000,sorted=FALSE)
  dmo<-finddmmm(expectboot=expectboot,expecttrue=expectboot,x=dp2m,interval=9,fast=TRUE,batch=10000,sorted=FALSE)
  all<-c(expectbootl=dlmo[1],etml=dlmo[2],ctml=dlmo[3],mx1l=dlmo[4],quatileexpectbootl=dlmo[5],
         expectboot=dmo[1],etm=dmo[2],ctm=dmo[3],mx1=dmo[4],quatileexpectboot=dmo[5]
  )
  return(all)
}

findd<-function(expectbootl=NULL,etml=NULL,ctml=NULL,mx1l=NULL,quatileexpectbootl=NULL,
                expectboot=NULL,etm=NULL,ctm=NULL,mx1=NULL,quatileexpectboot=NULL){
  quatileexpectbootl<-quatileexpectbootl
  etm1l<-c(1,etml,ctml)
  mx1l<-mx1l
  mx2l<-1/2
  if (mx1l>0.5){
    dqml1<-log(((quatileexpectbootl-mx1l)/(abs(1-mx1l)*((mx1l-mx2l)*2))),base=((abs(mx1l-mx2l)*2)))
  }else{
    quatileexpectbootl<-1-quatileexpectbootl
    mx1l<-1-mx1l
    dqml1<-log(((quatileexpectbootl-mx1l)/(abs(1-mx1l)*((mx1l-mx2l)*2))),base=((abs(mx1l-mx2l)*2)))
  }
  drml1<-(expectbootl-etm1l[2])/(etm1l[2]-etm1l[3])
  
  quatileexpectboot<-quatileexpectboot
  etm1<-c(1,etm,ctm)
  mx1<-mx1
  mx2<-1/2
  if (mx1>0.5){
    dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
  }else{
    quatileexpectboot<-1-quatileexpectboot
    quatileexpecttrue<-1-quatileexpecttrue
    mx1<-1-mx1
    dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
  }
  drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
  listd<-c(drml1,dqml1,drm1,dqm1)
  return(listd)
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

bootsize=54000
batchsizebase<-5400

#unbiased finite sample d is more difficult to process, first process each sample size, and then combine all together. Also, it is computationally very expensive to ensure a high accuracy.

library(Lmoments)
for(i in c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080))){
  batchsize<-round(batchsizebase/i)+1
  dmmm1boot<-c()
  dscale1boot<-c()
  dtm1boot<-c()
  dfm1boot<-c()
  for (j in (1:batchsize)){
  x<-c(rexp(i,1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  momentsx<-moments(x)
  lmomentsx<-Lmoments(x)
  dmmm<-finddmmm(expectboot=1,expecttrue=1,x,interval=9,fast=TRUE,batch=1000,sorted=TRUE)
  dmmm1boot<-rbind(dmmm1boot,dmmm)
  dscale1boot1<-finddscale(x,expectbootl=1/2,expectboot=1,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  dscale1boot1<-c(dscale1boot1,lmomentsx[2],momentsx[2],1/2,1)
  dscale1boot<-rbind(dscale1boot,dscale1boot1)
  dtm1boot1<-finddtm(x,expectbootl=1/6,expectboot=2,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  dtm1boot1<-c(dtm1boot1,lmomentsx[3],momentsx[3],1/6,2)
  dtm1boot<-rbind(dtm1boot,dtm1boot1)
  dfm1boot1<-finddfm(x,expectbootl=1/12,expectboot=9,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  dfm1boot1<-c(dfm1boot1,lmomentsx[4],momentsx[4],1/12,9)
  dfm1boot<-rbind(dfm1boot,dfm1boot1)}
  
  write.csv(dmmm1boot,paste("simulateddmmm_exp",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  colMeans(dscale1boot)
  write.csv(dscale1boot,paste("simulateddscale_exp",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dscale1boot<-colMeans(dscale1boot)
  
  findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
        expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
  dscale_exp<-findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
                    expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
  write.csv(dscale_exp,paste("dscale_exp",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  colMeans(dtm1boot)
  write.csv(dtm1boot,paste("simulateddtm_exp",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dtm1boot<-colMeans(dtm1boot)
  
  findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
        expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])
  
  dtm_exp<-findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
                 expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])
  
  write.csv(dtm_exp,paste("dtm_exp",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  colMeans(dfm1boot)
  write.csv(dfm1boot,paste("simulateddfm_exp",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dfm1boot<-colMeans(dfm1boot)
  
  findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
        expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])
  
  dfm_exp<-findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
                 expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])
  
  write.csv(dfm_exp,paste("dfm_exp",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
}

for(i in c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080))){
  batchsize<-round(batchsizebase/i)+1
  dmmm1boot<-c()
  dscale1boot<-c()
  dtm1boot<-c()
  dfm1boot<-c()
  for (j in (1:batchsize)){
    x<-rRayleigh(n=i,scale =1)
    x<-sort(x,decreasing = FALSE,method ="radix")
    momentsx<-moments(x)
    lmomentsx<-Lmoments(x)
    dmmm<-finddmmm(expectboot=1,expecttrue=1,x,interval=9,fast=TRUE,batch=1000,sorted=TRUE)
    dmmm1boot<-rbind(dmmm1boot,dmmm)
    dscale1boot1<-finddscale(x,expectbootl=0.5*(sqrt(2)-1)*sqrt(pi),expectboot=(2-(pi/2)),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
    dscale1boot1<-c(dscale1boot1,lmomentsx[2],momentsx[2],0.5*(sqrt(2)-1)*sqrt(pi),(2-(pi/2)))
    dscale1boot<-rbind(dscale1boot,dscale1boot1)
    dtm1boot1<-finddtm(x,expectbootl=(1/6)*(2*sqrt(6)+3*sqrt(2)-9)*sqrt(pi),expectboot=((sqrt(2-(pi/2)))^3)*2*sqrt((pi))*(pi-3)/((4-pi)^(3/2)),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
    dtm1boot1<-c(dtm1boot1,lmomentsx[3],momentsx[3],(1/6)*(2*sqrt(6)+3*sqrt(2)-9)*sqrt(pi),((sqrt(2-(pi/2)))^3)*2*sqrt((pi))*(pi-3)/((4-pi)^(3/2)))
    dtm1boot<-rbind(dtm1boot,dtm1boot1)
    dfm1boot1<-finddfm(x,expectbootl=(sqrt((77/6)-5*sqrt(6))-3/4)*sqrt(2*pi),expectboot=((sqrt(2-(pi/2)))^4)*(3-(6*(pi)^2-24*(pi)+16)/((4-pi)^(2))),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
    dfm1boot1<-c(dfm1boot1,lmomentsx[4],momentsx[4],(sqrt((77/6)-5*sqrt(6))-3/4)*sqrt(2*pi),((sqrt(2-(pi/2)))^4)*(3-(6*(pi)^2-24*(pi)+16)/((4-pi)^(2))))
    dfm1boot<-rbind(dfm1boot,dfm1boot1)}
  
  write.csv(dmmm1boot,paste("simulateddmmm_Rayleigh",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  colMeans(dscale1boot)
  write.csv(dscale1boot,paste("simulateddscale_Rayleigh",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dscale1boot<-colMeans(dscale1boot)
  
  findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
        expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
  dscale_Rayleigh<-findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
                         expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
  
  write.csv(dscale_Rayleigh,paste("dscale_Rayleigh",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  
  colMeans(dtm1boot)
  write.csv(dtm1boot,paste("simulateddtm_Rayleigh",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dtm1boot<-colMeans(dtm1boot)
  
  findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
        expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])
  
  dtm_Rayleigh<-findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
                      expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])
  
  write.csv(dtm_Rayleigh,paste("dtm_Rayleigh",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  colMeans(dfm1boot)
  write.csv(dfm1boot,paste("simulateddfm_Rayleigh",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dfm1boot<-colMeans(dfm1boot)
  
  findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
        expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])
  
  dfm_Rayleigh<-findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
                      expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])
  
  write.csv(dfm_Rayleigh,paste("dfm_Rayleigh",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  
}

for(i in c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080))){
  batchsize<-round(batchsizebase/i)+1
  dscale1boot<-c()
  dtm1boot<-c()
  dfm1boot<-c()
  for (j in (1:batchsize)){
    x<-rLaplace(n=i,location =0,scale=1)
    x<-sort(x,decreasing = FALSE,method ="radix")
    momentsx<-moments(x)
    lmomentsx<-Lmoments(x)
    dscale1boot1<-finddscale(x,expectbootl=3/4,expectboot=(2),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
    dscale1boot1<-c(dscale1boot1,lmomentsx[2],momentsx[2],3/4,(2))
    dscale1boot<-rbind(dscale1boot,dscale1boot1)
    dtm1boot1<-finddtm(x,expectbootl=0,expectboot=0,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
    dtm1boot1<-c(dtm1boot1,lmomentsx[3],momentsx[3],0,0)
    dtm1boot<-rbind(dtm1boot,dtm1boot1)
    dfm1boot1<-finddfm(x,expectbootl=(3/4)*1/(3*sqrt(2)),expectboot=6*(sqrt(2)^4),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
    dfm1boot1<-c(dfm1boot1,lmomentsx[4],momentsx[4],(3/4)*1/(3*sqrt(2)),6*(sqrt(2)^4))
    dfm1boot<-rbind(dfm1boot,dfm1boot1)}
  
  colMeans(dscale1boot)
  write.csv(dscale1boot,paste("simulateddscale_Laplace",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dscale1boot<-colMeans(dscale1boot)
  
  findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
        expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
  
  dscale_Laplace<-findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
                        expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
  
  write.csv(dscale_Laplace,paste("dscale_Laplace",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  colMeans(dtm1boot)
  write.csv(dtm1boot,paste("simulateddtm_Laplace",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dtm1boot<-colMeans(dtm1boot)
  
  findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
        expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])
  
  dtm_Laplace<-findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
                      expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])
  
  write.csv(dtm_Laplace,paste("dtm_Laplace",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  colMeans(dfm1boot)
  write.csv(dfm1boot,paste("simulateddfm_Laplace",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dfm1boot<-colMeans(dfm1boot)
  
  findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
        expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])
  
  dfm_Laplace<-findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
                     expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])
  
  write.csv(dfm_Laplace,paste("dfm_Laplace",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
}


for(i in c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080))){
  batchsize<-round(batchsizebase/i)+1
  dscale1boot<-c()
  dtm1boot<-c()
  dfm1boot<-c()
  for (j in (1:batchsize)){
    x<-rnorm(n=i,0,1)
    x<-sort(x,decreasing = FALSE,method ="radix")
    momentsx<-moments(x)
    lmomentsx<-Lmoments(x)
    dscale1boot1<-finddscale(x,expectbootl=1/sqrt(pi),expectboot=1,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
    dscale1boot1<-c(dscale1boot1,lmomentsx[2],momentsx[2],1/sqrt(pi),(1))
    dscale1boot<-rbind(dscale1boot,dscale1boot1)
    dtm1boot1<-finddtm(x,expectbootl=0,expectboot=0,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
    dtm1boot1<-c(dtm1boot1,lmomentsx[3],momentsx[3],0,0)
    dtm1boot<-rbind(dtm1boot,dtm1boot1)
    
    dfm1boot1<-finddfm(x,expectbootl=(1/sqrt(pi))*(30*(1/(pi))*(atan(sqrt(2)))-9),expectboot=3,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
    dfm1boot1<-c(dfm1boot1,lmomentsx[4],momentsx[4],(1/sqrt(pi))*(30*(1/(pi))*(atan(sqrt(2)))-9),3)
    dfm1boot<-rbind(dfm1boot,dfm1boot1)}
  
  colMeans(dscale1boot)
  write.csv(dscale1boot,paste("simulateddscale_norm",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dscale1boot<-colMeans(dscale1boot)
  
  findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
        expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
  
  dscale_norm<-findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
                     expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
  
  write.csv(dscale_norm,paste("dscale_norm",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  colMeans(dtm1boot)
  write.csv(dtm1boot,paste("simulateddtm_norm",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dtm1boot<-colMeans(dtm1boot)
  
  findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
        expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])
  
  dtm_norm<-findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
                     expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])
  
  write.csv(dtm_norm,paste("dtm_norm",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  
  colMeans(dfm1boot)
  write.csv(dfm1boot,paste("simulateddfm_norm",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dfm1boot<-colMeans(dfm1boot)
  
  findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
        expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])
  
  dfm_norm<-findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
                  expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])
  
  write.csv(dfm_norm,paste("dfm_norm",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
}


for(i in c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080))){
  batchsize<-round(batchsizebase/i)+1
  dscale1boot<-c()
  dtm1boot<-c()
  dfm1boot<-c()
  for (j in (1:batchsize)){
    x<-rlogis(n=i,location = 0, scale = 1)
    x<-sort(x,decreasing = FALSE,method ="radix")
    momentsx<-moments(x)
    lmomentsx<-Lmoments(x)
    dscale1boot1<-finddscale(x,expectbootl=1,expectboot=((pi^2)/3),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
    dscale1boot1<-c(dscale1boot1,lmomentsx[2],momentsx[2],1,(((pi^2)/3)))
    dscale1boot<-rbind(dscale1boot,dscale1boot1)
    dtm1boot1<-finddtm(x,expectbootl=0,expectboot=0,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
    dtm1boot1<-c(dtm1boot1,lmomentsx[3],momentsx[3],0,0)
    dtm1boot<-rbind(dtm1boot,dtm1boot1)
    
    dfm1boot1<-finddfm(x,expectbootl=1/6,expectboot=((6/5)+3)*((sqrt((pi^2)/3))^4),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
    dfm1boot1<-c(dfm1boot1,lmomentsx[4],momentsx[4],1/6,((6/5)+3)*((sqrt((pi^2)/3))^4))
    dfm1boot<-rbind(dfm1boot,dfm1boot1)}
  
  colMeans(dscale1boot)
  write.csv(dscale1boot,paste("simulateddscale_logis",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dscale1boot<-colMeans(dscale1boot)
  
  findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
        expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
  dscale_logis<-findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
                      expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
  
  write.csv(dscale_logis,paste("dscale_logis",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  colMeans(dtm1boot)
  write.csv(dtm1boot,paste("simulateddtm_logis",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dtm1boot<-colMeans(dtm1boot)
  
  findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
        expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])
  
  dtm_logis<-findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
                     expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])
  
  write.csv(dtm_logis,paste("dtm_logis",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
  colMeans(dfm1boot)
  write.csv(dfm1boot,paste("simulateddfm_logis",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  dfm1boot<-colMeans(dfm1boot)
  
  findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
        expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])
  
  dfm_logis<-findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
                   expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])
  
  write.csv(dfm_logis,paste("dfm_logis",batchnumber,i,".csv", sep = ","), row.names = TRUE)
  
}
