





scale<-function (x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE){
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
  all<-c(truel2=expectdps,bootl2=mean(dp2lm),truevar=expectdp2s,bootvar=mean(dp2m))
  return(all)
}
tm<-function (x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE){
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
  all<-c(truel3=expectdps,bootl3=mean(dp2lm),truetm=expectdp2s,boottm=mean(dp2m))
  return(all)
}

fm<-function (x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE){
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
  all<-c(truel4=expectdps,ebootl4=mean(dp2lm),truefm=expectdp2s,ebootfm=mean(dp2m))
  return(all)
}

#bootstrap works very good and the performance for the fourth moment is slighly better than sample kurtosis ~0.004, 
#because sample kurtosis is biased (finite sample bias), and uncorrectable without specific distribution assumption, 
#while U-statistics is much less biased.
library(Lmoments)

allfornorm<-c()
for(i in (1:10)){
  x<-c(rnorm(5400))
  x<-sort(x,decreasing = FALSE,method ="radix")
  scale1boot<-scale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  tm1<-tm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  fm1<-fm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  targetvar<-1
  targettm<-0
  targetfm<-3
  targetl2<-1/sqrt(pi)
  targetl3<-0
  targetl4<-(1/sqrt(pi))*(30*(1/(pi))*(atan(sqrt(2)))-9)
  all<-c(targetl2,targetvar,targetl3,targettm,targetl4,targetfm,scale1boot,tm1,fm1)
  
  allfornorm<-rbind(allfornorm,all)
}
allfornorm[is.infinite(allfornorm)] <-NA

write.csv(allfornorm,paste("Bootstrapnorm,",batchnumber,".csv", sep = ","), row.names = TRUE)

allforlaplace<-c()
for(i in (1:10)){
  library(VGAM)
  x<-c(rlaplace(5400, location = 0, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  scale1boot<-scale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  tm1<-tm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  fm1<-fm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  targetvar<-(2)
  targettm<-0
  targetfm<-6*(sqrt(2)^4)
  targetl2<-3/4
  targetl3<-0
  targetl4<-(3/4)*1/(3*sqrt(2))
  all<-c(targetl2,targetvar,targetl3,targettm,targetl4,targetfm,scale1boot,tm1,fm1)
  allforlaplace<-rbind(allforlaplace,all)
}
allforlaplace[is.infinite(allforlaplace)] <-NA

write.csv(allforlaplace,paste("Bootstraplaplace,",batchnumber,".csv", sep = ","), row.names = TRUE)


allforlogis<-c()
for(i in (1:10)){
  x<-c(rlogis(5400, location = 0, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  scale1boot<-scale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  tm1<-tm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  fm1<-fm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  targetvar<-((pi^2)/3)
  targettm<-0
  targetfm<-((6/5)+3)*((sqrt((pi^2)/3))^4)
  targetl2<-1
  targetl3<-0
  targetl4<-1/6
  all<-c(targetl2,targetvar,targetl3,targettm,targetl4,targetfm,scale1boot,tm1,fm1)
  
  allforlogis<-rbind(allforlogis,all)
}
allforlogis[is.infinite(allforlogis)] <-NA

write.csv(allforlogis,paste("Bootstraplogis,",batchnumber,".csv", sep = ","), row.names = TRUE)


allforRayleigh<-c()
for(i in (1:10)){
  library(VGAM)
  x<-c(rrayleigh(5400, scale = 1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  scale1boot<-scale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  tm1<-tm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  fm1<-fm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  targetvar<-(2-(pi/2))
  targettm<-((sqrt(2-(pi/2)))^3)*2*sqrt((pi))*(pi-3)/((4-pi)^(3/2))
  targetfm<-((sqrt(2-(pi/2)))^4)*(3-(6*(pi)^2-24*(pi)+16)/((4-pi)^(2)))
  targetl2<-0.5*(sqrt(2)-1)*sqrt(pi)
  targetl3<-(1/6)*(2*sqrt(6)+3*sqrt(2)-9)*sqrt(pi)
  targetl4<-(sqrt((77/6)-5*sqrt(6))-3/4)*sqrt(2*pi)
  all<-c(targetl2,targetvar,targetl3,targettm,targetl4,targetfm,scale1boot,tm1,fm1)
  
  allforRayleigh<-rbind(allforRayleigh,all)
}
allforRayleigh[is.infinite(allforRayleigh)] <-NA

write.csv(allforRayleigh,paste("BootstrapRayleigh,",batchnumber,".csv", sep = ","), row.names = TRUE)

allforexp<-c()
for(i in (1:10)){
  x<-c(rexp(5400,1))
  x<-sort(x,decreasing = FALSE,method ="radix")
  scale1boot<-scale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  tm1<-tm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  fm1<-fm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
  targetvar<-1
  targettm<-2
  targetfm<-9
  targetl2<-1/2
  targetl3<-1/6
  targetl4<-1/12
  all<-c(targetl2,targetvar,targetl3,targettm,targetl4,targetfm,scale1boot,tm1,fm1)
  allforexp<-rbind(allforexp,all)
}
allforexp[is.infinite(allforexp)] <-NA

write.csv(allforexp,paste("Bootstrapexp,",batchnumber,".csv", sep = ","), row.names = TRUE)

library(lmom)
listWeibull<-data.frame()
for (a in (1:300)) {
  allforweibull<-c()
  for(i in (1:10)){
    x<-c(rweibull(5400, shape=a/100, scale = 1))
    targetwei<-lmrwei(para = c(0, 1, a/100), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    scale1boot<-scale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    tm1<-tm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    fm1<-fm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    targetvar<-(gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2)
    targettm<-((sqrt(gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2))^3)*(gamma(1+3/(a/100))-3*(gamma(1+1/(a/100)))*((gamma(1+2/(a/100))))+2*((gamma(1+1/(a/100)))^3))/((sqrt(gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2))^(3))
    targetfm<-((sqrt(gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2))^4)*(gamma(1+4/(a/100))-4*(gamma(1+3/(a/100)))*((gamma(1+1/(a/100))))+6*(gamma(1+2/(a/100)))*((gamma(1+1/(a/100)))^2)-3*((gamma(1+1/(a/100)))^4))/(((gamma(1+2/(a/100))-(gamma(((1+1/(a/100)))))^2))^(2))
    targetl2<-targetwei[2]
    targetl3<-targetwei[3]*targetwei[2]
    targetl4<-targetwei[4]*targetwei[2]
    all<-c(targetl2,targetvar,targetl3,targettm,targetl4,targetfm,scale1boot,tm1,fm1)
    allforweibull<-rbind(allforweibull,all)
  }
  allforweibull[is.infinite(allforweibull)] <-NA
  listWeibull<-rbind(listWeibull,allforweibull)
}

write.csv(listWeibull,paste("BootstrapWeibull",batchnumber,".csv", sep = ","), row.names = TRUE)

library(lmom)
listgamma<-data.frame()
for (a in (1:300)) {
  allforgamma<-c()
  for(i in (1:10)){
    x<-c(rgamma(5400, shape=a/100, rate = 1))
    targetgam<-lmrgam(para = c(a/100, 1), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    scale1boot<-scale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    tm1<-tm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    fm1<-fm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    targetvar<-(a/100)
    targettm<-((sqrt(a/100))^3)*2/sqrt(a/100)
    targetfm<-((sqrt(a/100))^4)*((6/(a/100))+3)
    targetl2<-targetgam[2]
    targetl3<-targetgam[3]*targetgam[2]
    targetl4<-targetgam[4]*targetgam[2]
    all<-c(targetl2,targetvar,targetl3,targettm,targetl4,targetfm,scale1boot,tm1,fm1)
    allforgamma<-rbind(allforgamma,all)
  }
  allforgamma[is.infinite(allforgamma)] <-NA
  listgamma<-rbind(listgamma,allforgamma)
}

write.csv(listgamma,paste("Bootstrapgamma",batchnumber,".csv", sep = ","), row.names = TRUE)


listlnorm<-data.frame()
for (a in (1:300)) {
  allforlnorm<-c()
  for(i in (1:10)){
    x<-c(rlnorm(5400,meanlog=0,sdlog=a/100))
    targetlnorm<-lmrln3(para = c(0,0, a/100), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    scale1boot<-scale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    tm1<-tm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    fm1<-fm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    targetvar<-(exp((a/100)^2)*(-1+exp((a/100)^2)))
    targettm<-sqrt(exp((a/100)^2)-1)*((2+exp((a/100)^2)))*((sqrt(exp((a/100)^2)*(-1+exp((a/100)^2))))^3)
    targetfm<-((-3+exp(4*((a/100)^2))+2*exp(3*((a/100)^2))+3*exp(2*((a/100)^2))))*((sqrt(exp((a/100)^2)*(-1+exp((a/100)^2))))^4)
    targetl2<-targetlnorm[2]
    targetl3<-targetlnorm[3]*targetlnorm[2]
    targetl4<-targetlnorm[4]*targetlnorm[2]
    all<-c(targetl2,targetvar,targetl3,targettm,targetl4,targetfm,scale1boot,tm1,fm1)
    allforlnorm<-rbind(allforlnorm,all)
  }
  allforlnorm[is.infinite(allforlnorm)] <-NA
  listlnorm<-rbind(listlnorm,allforlnorm)
}

write.csv(listlnorm,paste("Bootstraplnorm",batchnumber,".csv", sep = ","), row.names = TRUE)
#notice that third moment only defined for shape>3, fourth moment only defined for shape>4
listpareto<-data.frame()
for (a in (1:300)) {
  allforpareto<-c()
  for(i in (1:10)){
    library(VGAM)
    x<-c(VGAM::rpareto(5400, scale  = 1, shape=2+a/100))
    targetlpareto<-lmrgpa(para = c(1,1/(2+a/100),- 1/(2+a/100)), nmom = 4)
    x<-sort(x,decreasing = FALSE,method ="radix")
    scale1boot<-scale(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    tm1<-tm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    fm1<-fm(x,interval=9,fast=TRUE,batch=1000,boot=TRUE,subsample=54000,sorted=TRUE)
    targetvar<-(((2+a/100))*(1)/((-2+(2+a/100))*((-1+(2+a/100))^2)))
    targettm<-((((2+a/100)+1)*(2)*(sqrt(a/100)))/((-3+(2+a/100))*(((2+a/100))^(1/2))))*(((sqrt(((2+a/100))*(1)/((-2+(2+a/100))*((-1+(2+a/100))^2))))^3))
    targetfm<-(3+(6*((2+a/100)^3+(2+a/100)^2-6*(2+a/100)-2)/(((2+a/100))*((-3+(2+a/100)))*((-4+(2+a/100))))))*((sqrt(((2+a/100))*(1)/((-2+(2+a/100))*((-1+(2+a/100))^2))))^4)
    targetl2<-targetlpareto[2]
    targetl3<-targetlpareto[3]*targetlpareto[2]
    targetl4<-targetlpareto[4]*targetlpareto[2]
    all<-c(targetl2,targetvar,targetl3,targettm,targetl4,targetfm,scale1boot,tm1,fm1)
    
    allforpareto<-rbind(allforpareto,all)
  }
  allforpareto[is.infinite(allforpareto)] <-NA
  listpareto<-rbind(listpareto,allforpareto)
}

write.csv(listpareto,paste("Bootstrappareto",batchnumber,".csv", sep = ","), row.names = TRUE)
