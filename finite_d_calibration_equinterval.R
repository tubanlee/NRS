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
  quatileexpectboot<-(min(which(sortedx>(expectboot)))-1)/length(x)
  quatileexpecttrue<-(min(which(sortedx>(expecttrue)))-1)/length(x)
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
  lm1<-Lmoments(sortedx)
  
  expectdps<-lm1[3]
  expectdp2s<-((sum((sortedx - mean(sortedx))^3)/lengthn)*(lengthn^2/((lengthn-1)*(lengthn-2))))
  dlmo<-finddmmm(expectboot=expectbootl,expecttrue=expectbootl,x=dp2lm,interval=9,fast=fast,batch=batch,sorted=FALSE)
  dmo<-finddmmm(expectboot=expectboot,expecttrue=expectboot,x=dp2m,interval=9,fast=fast,batch=batch,sorted=FALSE)
  all<-c(drl3=dlmo[1],dql3=dlmo[2],
         drtm=dmo[1],dqtm=dmo[2]
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
  lm1<-Lmoments(sortedx)
  expectdps<-lm1[4]
  expectdp2s<-(sum((sortedx - mean(sortedx))^4)/lengthn)
  dlmo<-finddmmm(expectboot=expectbootl,expecttrue=expectbootl,x=dp2lm,interval=9,fast=fast,batch=batch,sorted=FALSE)
  dmo<-finddmmm(expectboot=expectboot,expecttrue=expectboot,x=dp2m,interval=9,fast=fast,batch=batch,sorted=FALSE)
  all<-c(drl4=dlmo[1],dql4=dlmo[2],
         drfm=dmo[1],dqfm=dmo[2]  )
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
  lm1<-Lmoments(sortedx)
  expectdps<-lm1[2]
  expectdp2s<-(sd(sortedx))^2
  dlmo<-finddmmm(expectboot=expectbootl,expecttrue=expectbootl,x=dp2lm,interval=9,fast=TRUE,batch=10000,sorted=FALSE)
  dmo<-finddmmm(expectboot=expectboot,expecttrue=expectboot,x=dp2m,interval=9,fast=TRUE,batch=10000,sorted=FALSE)
  all<-c(drl2=dlmo[1],dql2=dlmo[2],
         drvar=dmo[1],dqvar=dmo[2]   )
  return(all)
}

eRayleigh<-function (n, scale) {
  sample1 <- scale * sqrt(-2 * log((seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))))
  sample1[scale <= 0] <- NaN
  sample1
}
eexp<-function (n, scale) {
  sample1 <- (-1/scale)*(log(1-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))))
  sample1[scale <= 0] <- NaN
  sample1
}
enorm<-function (n, location,scale) {
  library(pracma)
  sample1 <- location+(scale)*sqrt(2)*erfinv(2*(1-seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))-1)
  sample1
}
eLaplace<-function (n,location,scale) {
  sample1<-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))
  sample1<-location - sign(sample1 - 0.5) * scale * (log(2) + ifelse(sample1 < 0.5, log(sample1), log1p(-sample1)))
  sample1
}
elogis<-function (n,location,scale) {
  sample1<-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))
  sample1<-location + scale * log((1-sample1)/sample1)
  sample1
}
library(foreach)
library(doParallel)
numCores <- detectCores()
registerDoParallel(numCores) 
#parallel
bootsize=216000
simulatedbatchLaplace<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
  library(Lmoments)
  x<-c(eLaplace(n=i, location = 0, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  dmmm<-c(0,1000000)
  dscale1boot<-finddscale(x,expectbootl=3/4,expectboot=2,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  dtm1boot<-c(0,1000000,0,1000000)
  dfm1boot<-finddfm(x,expectbootl=(3/4)*1/(3*sqrt(2)),expectboot=6*(sqrt(2)^4),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  all<-c(dmmm,dscale1boot,dtm1boot,dfm1boot)
}

simulatedbatchLaplace[is.infinite(simulatedbatchLaplace)] <-NA
write.csv(simulatedbatchLaplace,paste("fdequinterval_Laplace.csv", sep = ","), row.names = TRUE)


simulatedbatchnorm<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
  library(Lmoments)
  x<-c(enorm(n=i,location=0,scale=1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  dmmm<-c(0,1000000)
  dscale1boot<-finddscale(x,expectbootl=1/sqrt(pi),expectboot=1,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  dtm1boot<-c(0,1000000,0,1000000)
  dfm1boot<-finddfm(x,expectbootl=(1/sqrt(pi))*(30*(1/(pi))*(atan(sqrt(2)))-9),expectboot=3,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  all<-c(dmmm,dscale1boot,dtm1boot,dfm1boot)
}

simulatedbatchnorm[is.infinite(simulatedbatchnorm)] <-NA
write.csv(simulatedbatchnorm,paste("fdequinterval_norm.csv", sep = ","), row.names = TRUE)

simulatedbatchlogis<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
  library(Lmoments)
  x<-c(elogis(n=i, location = 0, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  dmmm<-c(0,1000000)
  dscale1boot<-finddscale(x,expectbootl=1,expectboot=((pi^2)/3),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  dtm1boot<-c(0,1000000,0,1000000)
  dfm1boot<-finddfm(x,expectbootl=1/6,expectboot=((6/5)+3)*((sqrt((pi^2)/3))^4),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  all<-c(dmmm,dscale1boot,dtm1boot,dfm1boot)
}

simulatedbatchlogis[is.infinite(simulatedbatchlogis)] <-NA
write.csv(simulatedbatchlogis,paste("fdequinterval_logis.csv", sep = ","), row.names = TRUE)

simulatedbatch<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
  library(Lmoments)
  x<-c(eexp(n=i,scale=1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  dmmm<-finddmmm(expectboot=1,expecttrue=1,x,interval=9,fast=TRUE,batch=1000,sorted=TRUE)
  dscale1boot<-finddscale(x,expectbootl=1/2,expectboot=1,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  dtm1boot<-finddtm(x,expectbootl=1/6,expectboot=2,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  dfm1boot<-finddfm(x,expectbootl=1/12,expectboot=9,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  all<-c(dmmm,dscale1boot,dtm1boot,dfm1boot)
}

simulatedbatch[is.infinite(simulatedbatch)] <-NA
write.csv(simulatedbatch,paste("fdequinterval_exp.csv", sep = ","), row.names = TRUE)

simulatedbatchRayleigh<-foreach(i = c(seq(from=9, to=108, by=1),seq(from=117, to=1305, by=108),seq(from=1413, to=5733, by=1080)), .combine = 'rbind') %dopar% {
  library(Lmoments)
  x<-c(eRayleigh(n=i,scale=1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  dmmm<-finddmmm(expectboot=sqrt(pi/2),expecttrue=sqrt(pi/2),x,interval=9,fast=TRUE,batch=1000,sorted=TRUE)
  dscale1boot<-finddscale(x,expectbootl=0.5*(sqrt(2)-1)*sqrt(pi),expectboot=(2-(pi/2)),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  dtm1boot<-finddtm(x,expectbootl=(1/6)*(2*sqrt(6)+3*sqrt(2)-9)*sqrt(pi),expectboot=((sqrt(2-(pi/2)))^3)*2*sqrt((pi))*(pi-3)/((4-pi)^(3/2)),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  dfm1boot<-finddfm(x,expectbootl=(sqrt((77/6)-5*sqrt(6))-3/4)*sqrt(2*pi),expectboot=((sqrt(2-(pi/2)))^4)*(3-(6*(pi)^2-24*(pi)+16)/((4-pi)^(2))),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=TRUE)
  all<-c(dmmm,dscale1boot,dtm1boot,dfm1boot)
}

simulatedbatchRayleigh[is.infinite(simulatedbatchRayleigh)] <-NA
write.csv(simulatedbatchRayleigh,paste("fdequinterval_Rayleigh.csv", sep = ","), row.names = TRUE)

