


#the goal of this file is to estimate the asymptotic d for five single parameter distributions with accuracy to three decimal places, i.e..
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
    groupp<-foreach (i=0:8, .combine=c) %dopar% {
      sum(x_ordered[(i*IntKsamples+1):((i+1)*(IntKsamples))])
    }
    etmsum<-(sum(groupp[2],groupp[5],groupp[8]))
    ctmsum<-(sum(groupp[3],groupp[4],groupp[6],groupp[7]))
    Groupsumweight<-c(sum(etmsum,ctmsum,groupp[9],groupp[1]),sum(etmsum,groupp[2],groupp[8]),ctmsum)
    etmlength<-IntKsamples*5
    ctmlength<-IntKsamples*4
    Groupmean<-c(Groupsumweight[1]/(IntKsamples*9),Groupsumweight[2]/etmlength,Groupsumweight[3]/ctmlength)
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
    quatileexpecttrue<-1-quatileexpecttrue
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
ePareto<-function (n, scale, shape) {
  sample1 <- scale*((seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1))))^(-1/shape)
  sample1[scale <= 0] <- NaN
  sample1[shape <= 0] <- NaN
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
elnorm<-function (n,location,scale) {
  library(pracma)
  sample1 <- exp(location+sqrt(scale^2)*sqrt(2)*erfinv(2*(1-seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))-1))
  sample1
}
egamma<-function (n,shape,scale = 1) {
  sample1<-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))
  sample1<-qgamma(sample1,shape=shape,scale=scale)
  sample1
}
eWeibull<-function (n,shape, scale = 1){
  sample1<-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))
  sample1<-qweibull(sample1,shape=shape, scale = scale)
  sample1
}
factdivide<-function(n1,n2){
  decin1<-n1-floor(n1)
  if(decin1==0){decin1=1}
  decin2<-n2-floor(n2)
  if(decin2==0){decin2=1}
  n1seq<-seq(decin1,n1, by=1)
  n2seq<-seq(decin2,n2,by=1)
  all<-list(n1seq,n2seq)
  maxlen <- max(lengths(all))
  all2 <- as.data.frame(lapply(all, function(lst) c(lst, rep(1, maxlen - length(lst)))))
  division<-all2[,1]/all2[,2]
  answer<-exp(sum(log(division)))*(gamma(decin1)/gamma(decin2))
  return(answer)
}


x<-eRayleigh(n=9000,scale =1)

mean(x)
#the bias of equinterval simulation
abs(1.253193-sqrt(pi/2))
sd(x)
n=9000
#an approximation formula to estimate the bias
sd(x)*(1-(sqrt(2/(n-1)))*factdivide(n1=(n/2)-3.5,n2=((n-1)/2)-3.5))

x<-eRayleigh(n=90000,scale =1)

mean(x)
abs(1.253299-sqrt(pi/2))
sd(x)
n=90000
sd(x)*(1-(sqrt(2/(n-1)))*factdivide(n1=(n/2)-3.5,n2=((n-1)/2)-3.5))

x<-eRayleigh(n=900000,scale =1)

mean(x)
abs(1.253312-sqrt(pi/2))
sd(x)
n=900000
sd(x)*(1-(sqrt(2/(n-1)))*factdivide(n1=(n/2)-3.5,n2=((n-1)/2)-3.5))

x<-eRayleigh(n=9000000,scale =1)

mean(x)
abs(1.253314-sqrt(pi/2))
sd(x)
n=9000000
sd(x)*(1-(sqrt(2/(n-1)))*factdivide(n1=(n/2)-3.5,n2=((n-1)/2)-3.5))


library(Lmoments)

x<-eRayleigh(n=9000,scale =1)
Lmoments(x)
moments(x)
#the bias of equinterval simulation
#var
abs(0.4285605-((2-(pi/2))))
#tm
abs(0.1755594 -((sqrt(2-(pi/2)))^3)*2*sqrt((pi))*(pi-3)/((4-pi)^(3/2)))
#fm
abs(0.5900809 -((sqrt(2-(pi/2)))^4)*(3-(6*(pi)^2-24*(pi)+16)/((4-pi)^(2))))

#l2
abs(0.3669501-0.5*(sqrt(2)-1)*sqrt(pi))
#l3
abs(0.04174286 -(1/6)*(2*sqrt(6)+3*sqrt(2)-9)*sqrt(pi))
#l4
abs(0.03858559 -(sqrt((77/6)-5*sqrt(6))-3/4)*sqrt(2*pi))

x<-eRayleigh(n=90000,scale =1)
#the bias of equinterval simulation
#var
abs(0.4291203 -((2-(pi/2))))
#tm
abs(0.1771689  -((sqrt(2-(pi/2)))^3)*2*sqrt((pi))*(pi-3)/((4-pi)^(3/2)))
#fm
abs(0.5965140  -((sqrt(2-(pi/2)))^4)*(3-(6*(pi)^2-24*(pi)+16)/((4-pi)^(2))))

#l2
abs(0.3670709 -0.5*(sqrt(2)-1)*sqrt(pi))
#l3
abs(0.04182399  -(1/6)*(2*sqrt(6)+3*sqrt(2)-9)*sqrt(pi))
#l4
abs(0.03866775 -(sqrt((77/6)-5*sqrt(6))-3/4)*sqrt(2*pi))

x<-eRayleigh(n=9000000,scale =1)
Lmoments(x)
moments(x)
#the bias of equinterval simulation
#var
abs(0.4292024 -((2-(pi/2))))
#tm
abs(0.1774547  -((sqrt(2-(pi/2)))^3)*2*sqrt((pi))*(pi-3)/((4-pi)^(3/2)))
#fm
abs(0.5977691  -((sqrt(2-(pi/2)))^4)*(3-(6*(pi)^2-24*(pi)+16)/((4-pi)^(2))))

#l2
abs(0.367087 -0.5*(sqrt(2)-1)*sqrt(pi))
#l3
abs(0.04183571  -(1/6)*(2*sqrt(6)+3*sqrt(2)-9)*sqrt(pi))
#l4
abs(0.03867962 -(sqrt((77/6)-5*sqrt(6))-3/4)*sqrt(2*pi))

#so, with sample size 9 million is enough to ensure three decimal accuracy.


#require foreach and doparallel for parallel processing of bootstrap (not available for some types of computers)
if (!require("foreach")) install.packages("foreach")
library(foreach)
if (!require("doSNOW")) install.packages("doSNOW")
library(doSNOW)
if (!require("doParallel")) install.packages("doParallel")
library(doParallel)
#registering clusters, can set a smaller number using numCores-1 
numCores <- detectCores()
cl <- makeCluster(numCores, type = "SOCK")
registerDoSNOW(cl)
x<-c(eexp(n=9000000,scale=1))
mean(x)
abs(0.9999991-1)

dmmm<-finddmmme(expectboot=mean(x),expecttrue=mean(x),x,interval=9,fast=TRUE,batch=1000,sorted=TRUE)
dmmm
#the result is 0.3669145 and 0.8215032 , four decimal accuracy
#the aymptotic value of robust mean
0.366919 
abs(0.366919-0.3669145)
#quantile mean
0.821497
abs(0.821497-0.8215032)

#further increase the sample size to 21600000 (the memory limit of most PC)

#the analytical solution is shown in SI, which is very complex and time consuming, also, for some especial distributions, the analytical solution is quite remote.
x<-c(eexp(n=21600000,scale=1))
mean(x)
abs(0.9999996-1)

dmmm<-finddmmme(expectboot=mean(x),expecttrue=mean(x),x,interval=9,fast=TRUE,batch=1000,sorted=TRUE)
dmmm
#the result is 0.3669145 and 0.8215032 , five decimal accuracy
#the aymptotic value of robust mean
0.366919 
abs(0.366919-0.3669170)
#quantile mean
0.821497
abs(0.821497-0.8215006)

samplesize=21600000
bootsize=21600
batchsize=20
dscale1boot<-c()
for(i in 1:batchsize){
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
      quatileexpecttrue<-1-quatileexpecttrue
      mx1<-1-mx1
      dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
    }
    drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
    listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
    return(listd)
  }
  x<-c(rexp(n=samplesize,1))
  dscale1boot1<-finddscale(x,expectbootl=1/2,expectboot=1,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=FALSE)
  dscale1boot<-rbind(dscale1boot,dscale1boot1)
}

colMeans(dscale1boot)
write.csv(dscale1boot,paste("simulateddscale_exp.csv", sep = ","), row.names = TRUE)
dscale1boot<-colMeans(dscale1boot)

findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
                expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
dscale_exp<-findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
                  expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])

write.csv(dscale_exp,paste("dscale_exp.csv", sep = ","), row.names = TRUE)


dtm1boot<-c()
for(i in 1:batchsize){
  eexp<-function (n, scale) {
    sample1 <- (-1/scale)*(log(1-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))))
    sample1[scale <= 0] <- NaN
    sample1
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
      quatileexpecttrue<-1-quatileexpecttrue
      mx1<-1-mx1
      dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
    }
    drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
    listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
    return(listd)
  }
  
  x<-c(rexp(n=samplesize,1))
  dtm1boot1<-finddtm(x,expectbootl=1/6,expectboot=2,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=FALSE)
  dtm1boot<-rbind(dtm1boot,dtm1boot1)
}

colMeans(dtm1boot)
write.csv(dtm1boot,paste("simulateddtm_exp.csv", sep = ","), row.names = TRUE)
dtm1boot<-colMeans(dtm1boot)

findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
      expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])

dtm_exp<-findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
               expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])

write.csv(dtm_exp,paste("dtm_exp.csv", sep = ","), row.names = TRUE)


dfm1boot<-c()
for(i in 1:batchsize){
  eexp<-function (n, scale) {
    sample1 <- (-1/scale)*(log(1-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))))
    sample1[scale <= 0] <- NaN
    sample1
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
      quatileexpecttrue<-1-quatileexpecttrue
      mx1<-1-mx1
      dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
    }
    drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
    listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
    return(listd)
  }
  x<-c(rexp(n=samplesize,1))
  dfm1boot1<-finddfm(x,expectbootl=1/12,expectboot=9,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=FALSE)
  dfm1boot<-rbind(dfm1boot,dfm1boot1)
}

colMeans(dfm1boot)
write.csv(dfm1boot,paste("simulateddfm_exp.csv", sep = ","), row.names = TRUE)
dfm1boot<-colMeans(dfm1boot)

findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
      expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])

dfm_exp<-findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
               expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])

write.csv(dfm_exp,paste("dfm_exp.csv", sep = ","), row.names = TRUE)


#the analytical solution is shown in SI, which is very complex and time consuming, also, for some especial distributions, the analytical solution is quite remote.
x<-eRayleigh(n=samplesize,scale =1)
mean(x)
abs(1.253314-sqrt(pi/2))

dmmm<-finddmmme(expectboot=mean(x),expecttrue=mean(x),x=x,interval=9,fast=TRUE,batch=1000,sorted=FALSE)
dmmm
#0.4029587 0.4446060
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

dscale1boot<-c()
for(i in 1:batchsize){
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
      quatileexpecttrue<-1-quatileexpecttrue
      mx1<-1-mx1
      dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
    }
    drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
    listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
    return(listd)
  }
  x<-c(rRayleigh(n=samplesize, scale = 1))
  dscale1boot1<-finddscale(x,expectbootl=0.5*(sqrt(2)-1)*sqrt(pi),expectboot=(2-(pi/2)),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=FALSE)
  dscale1boot<-rbind(dscale1boot,dscale1boot1)
}
colMeans(dscale1boot)
write.csv(dscale1boot,paste("simulateddscale_Rayleigh.csv", sep = ","), row.names = TRUE)
dscale1boot<-colMeans(dscale1boot)

findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
      expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
dscale_Rayleigh<-findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
               expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])

write.csv(dscale_Rayleigh,paste("dscale_Rayleigh.csv", sep = ","), row.names = TRUE)



dtm1boot<-c()
for(i in 1:batchsize){
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
      quatileexpecttrue<-1-quatileexpecttrue
      mx1<-1-mx1
      dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
    }
    drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
    listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
    return(listd)
  }
  
  x<-c(rRayleigh(n=samplesize, scale = 1))
  dtm1boot1<-finddtm(x,expectbootl=(1/6)*(2*sqrt(6)+3*sqrt(2)-9)*sqrt(pi),expectboot=((sqrt(2-(pi/2)))^3)*2*sqrt((pi))*(pi-3)/((4-pi)^(3/2)),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=FALSE)
  dtm1boot<-rbind(dtm1boot,dtm1boot1)
}

colMeans(dtm1boot)
write.csv(dtm1boot,paste("simulateddtm_Rayleigh.csv", sep = ","), row.names = TRUE)
dtm1boot<-colMeans(dtm1boot)

findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
      expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])

dtm_Rayleigh<-findd(expectbootl=dtm1boot[1],etml=dtm1boot[2],ctml=dtm1boot[3],mx1l=dtm1boot[4],quatileexpectbootl=dtm1boot[5],
                    expectboot=dtm1boot[6],etm=dtm1boot[7],ctm=dtm1boot[8],mx1=dtm1boot[9],quatileexpectboot=dtm1boot[10])

write.csv(dtm_Rayleigh,paste("dtm_Rayleigh.csv", sep = ","), row.names = TRUE)


dfm1boot<-c()
for(i in 1:batchsize){
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
      quatileexpecttrue<-1-quatileexpecttrue
      mx1<-1-mx1
      dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
    }
    drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
    listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
    return(listd)
  }
  x<-c(rRayleigh(n=samplesize, scale = 1))
  dfm1boot1<-finddfm(x,expectbootl=(sqrt((77/6)-5*sqrt(6))-3/4)*sqrt(2*pi),expectboot=((sqrt(2-(pi/2)))^4)*(3-(6*(pi)^2-24*(pi)+16)/((4-pi)^(2))),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=FALSE)
  dfm1boot<-rbind(dfm1boot,dfm1boot1)
}

colMeans(dfm1boot)
write.csv(dfm1boot,paste("simulateddfm_Rayleigh.csv", sep = ","), row.names = TRUE)
dfm1boot<-colMeans(dfm1boot)

findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
      expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])

dfm_Rayleigh<-findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
                    expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])

write.csv(dfm_Rayleigh,paste("dfm_Rayleigh.csv", sep = ","), row.names = TRUE)


x<-c(eLaplace(n=samplesize, location = 0, scale = 1))
mean(x)
abs(-1.273414e-15-0)
#dmmm
#0 0


dscale1boot<-c()
for(i in 1:batchsize){
  
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
      quatileexpecttrue<-1-quatileexpecttrue
      mx1<-1-mx1
      dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
    }
    drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
    listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
    return(listd)
  }
  x<-c(rLaplace(n=samplesize, location = 0, scale = 1))
  dscale1boot1<-finddscale(x,expectbootl=3/4,expectboot=(2),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=FALSE)
  dscale1boot<-rbind(dscale1boot,dscale1boot1)
}
colMeans(dscale1boot)
write.csv(dscale1boot,paste("simulateddscale_Laplace.csv", sep = ","), row.names = TRUE)
dscale1boot<-colMeans(dscale1boot)

findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
      expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])

dscale_Laplace<-findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
                      expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])

write.csv(dscale_Laplace,paste("dscale_Laplace.csv", sep = ","), row.names = TRUE)


#dtm1boot
#0 0

dfm1boot<-c()
for(i in 1:batchsize){
  
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
      quatileexpecttrue<-1-quatileexpecttrue
      mx1<-1-mx1
      dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
    }
    drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
    listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
    return(listd)
  }
  x<-c(rLaplace(n=samplesize, location = 0, scale = 1))
  dfm1boot1<-finddfm(x,expectbootl=(3/4)*1/(3*sqrt(2)),expectboot=6*(sqrt(2)^4),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=FALSE)
  dfm1boot<-rbind(dfm1boot,dfm1boot1)
}

colMeans(dfm1boot)
write.csv(dfm1boot,paste("simulateddfm_Laplace.csv", sep = ","), row.names = TRUE)
dfm1boot<-colMeans(dfm1boot)

findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
      expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])

dfm_Laplace<-findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
                   expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])

write.csv(dfm_Laplace,paste("dfm_Laplace.csv", sep = ","), row.names = TRUE)


x<-c(enorm(n=samplesize,location=0,scale=1))
mean(x)
abs(4.580853e-16-0)

#dmmm
#0 0


dscale1boot<-c()
for(i in 1:batchsize){
  
  finddscale<-function (x,expectbootl,expectboot,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE){
    if(sorted){
      sortedx<-x
    }else{
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
    }
    if (boot){
      subtract<-foreach (i=1:subsample, .combine=rbind) %dopar% {
        sort(sample(sortedx, size = 2))
      }
      dp2lm<-foreach (i=1:subsample, .combine=rbind) %dopar% {
        getlm<-function(vector){ 
          (vector[2]-vector[1])/2
        }
        getlm(subtract[i,])
      }
      dp2m<-foreach (i=1:subsample, .combine=rbind) %dopar% {
        getm<-function(vector){ 
          ((vector[1]-vector[2])^2)/2
        }
        getm(subtract[i,])
      }
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
      quatileexpecttrue<-1-quatileexpecttrue
      mx1<-1-mx1
      dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
    }
    drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
    listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
    return(listd)
  }
  x<-c(rnorm(n=samplesize,0,1))
  dscale1boot1<-finddscale(x,expectbootl=1/sqrt(pi),expectboot=1,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=FALSE)
  dscale1boot<-rbind(dscale1boot,dscale1boot1)
}

colMeans(dscale1boot)
write.csv(dscale1boot,paste("simulateddscale_norm.csv", sep = ","), row.names = TRUE)
dscale1boot<-colMeans(dscale1boot)

findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
      expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])

dscale_norm<-findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
                   expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])

write.csv(dscale_norm,paste("dscale_norm.csv", sep = ","), row.names = TRUE)


#dtm1boot
#0 0

dfm1boot<-c()
for(i in 1:batchsize){
  
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
      quatileexpecttrue<-1-quatileexpecttrue
      mx1<-1-mx1
      dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
    }
    drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
    listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
    return(listd)
  }
  x<-c(rnorm(n=samplesize,0,1))
  dfm1boot1<-finddfm(x,expectbootl=(1/sqrt(pi))*(30*(1/(pi))*(atan(sqrt(2)))-9),expectboot=3,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=FALSE)
  dfm1boot<-rbind(dfm1boot,dfm1boot1)
}

colMeans(dfm1boot)
write.csv(dfm1boot,paste("simulateddfm_norm.csv", sep = ","), row.names = TRUE)
dfm1boot<-colMeans(dfm1boot)

findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
      expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])

dfm_norm<-findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
                expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])

write.csv(dfm_norm,paste("dfm_norm.csv", sep = ","), row.names = TRUE)


x<-c(elogis(n=samplesize, location = 0, scale = 1))
mean(x)
abs(1.395458e-15-0)

#dmmm
#0 0


dscale1boot<-c()
for(i in 1:batchsize){
  
  finddscale<-function (x,expectbootl,expectboot,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE){
    if(sorted){
      sortedx<-x
    }else{
      sortedx<-sort(x,decreasing = FALSE,method ="radix")
    }
    if (boot){
      subtract<-foreach (i=1:subsample, .combine=rbind) %dopar% {
        sort(sample(sortedx, size = 2))
      }
      dp2lm<-foreach (i=1:subsample, .combine=rbind) %dopar% {
        getlm<-function(vector){ 
          (vector[2]-vector[1])/2
        }
        getlm(subtract[i,])
      }
      dp2m<-foreach (i=1:subsample, .combine=rbind) %dopar% {
        getm<-function(vector){ 
          ((vector[1]-vector[2])^2)/2
        }
        getm(subtract[i,])
      }
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
      quatileexpecttrue<-1-quatileexpecttrue
      mx1<-1-mx1
      dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
    }
    drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
    listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
    return(listd)
  }
  x<-c(rlogis(n=samplesize, location = 0, scale = 1))
  dscale1boot1<-finddscale(x,expectbootl=1,expectboot=((pi^2)/3),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=FALSE)
  dscale1boot<-rbind(dscale1boot,dscale1boot1)
}

colMeans(dscale1boot)
write.csv(dscale1boot,paste("simulateddscale_logis.csv", sep = ","), row.names = TRUE)
dscale1boot<-colMeans(dscale1boot)

findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
      expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])
dscale_logis<-findd(expectbootl=dscale1boot[1],etml=dscale1boot[2],ctml=dscale1boot[3],mx1l=dscale1boot[4],quatileexpectbootl=dscale1boot[5],
                expectboot=dscale1boot[6],etm=dscale1boot[7],ctm=dscale1boot[8],mx1=dscale1boot[9],quatileexpectboot=dscale1boot[10])

write.csv(dscale_logis,paste("dscale_logis.csv", sep = ","), row.names = TRUE)


#dtm1boot
#0 0
dfm1boot<-c()
for(i in 1:batchsize){
  elogis<-function (n,location,scale) {
    sample1<-(seq(from=1/(n+1), to=1-1/(n+1), by=1/(n+1)))
    sample1<-location + scale * log((1-sample1)/sample1)
    sample1
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
      quatileexpecttrue<-1-quatileexpecttrue
      mx1<-1-mx1
      dqm1<-log(((quatileexpectboot-mx1)/(abs(1-mx1)*((mx1-mx2)*2))),base=((abs(mx1-mx2)*2)))
    }
    drm1<-(expectboot-etm1[2])/(etm1[2]-etm1[3])
    listd<-c(expectboot,etm1[2],etm1[3],mx1,quatileexpectboot)
    return(listd)
  }
  x<-c(rlogis(n=samplesize, location = 0, scale = 1))
  dfm1boot1<-finddfm(x,expectbootl=1/6,expectboot=((6/5)+3)*((sqrt((pi^2)/3))^4),interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=bootsize,sorted=FALSE)
  dfm1boot<-rbind(dfm1boot,dfm1boot1)
}

colMeans(dfm1boot)
write.csv(dfm1boot,paste("simulateddfm_logis.csv", sep = ","), row.names = TRUE)
dfm1boot<-colMeans(dfm1boot)

findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
      expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])

dfm_logis<-findd(expectbootl=dfm1boot[1],etml=dfm1boot[2],ctml=dfm1boot[3],mx1l=dfm1boot[4],quatileexpectbootl=dfm1boot[5],
                 expectboot=dfm1boot[6],etm=dfm1boot[7],ctm=dfm1boot[8],mx1=dfm1boot[9],quatileexpectboot=dfm1boot[10])

write.csv(dfm_logis,paste("dfm_logis.csv", sep = ","), row.names = TRUE)
