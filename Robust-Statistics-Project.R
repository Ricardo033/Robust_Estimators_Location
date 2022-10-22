#============================================#
# Robust Statistics                          #
# MAD, Sn, Qn location estimators - projec 1 #
# Ricardo Castañeda - r0731529               #
# QingYuan Xue - r0829510                    #
# AL-AWA Sarah - r0688008                    #
# Zachary Jones - r0772334                   #                                          #
#============================================#

## Packages:

rm(list=ls())

#install.packages("DescTools")
#install.packages("reshape")
library(robustbase)
library(reshape)
library(ggplot2)
library(DescTools)
library(ggpubr)
suppressMessages(library(matrixStats))
library(reshape2)
library(reshape2, data.table)
library(dplyr)
library(latex2exp)
library(covSde)

#==============#
# Question 2   #
#==============#

## Functions
### function for Bickel-Lehmann estimator
BL<-function(x){
  len<-length(x)
  lenn<-len*(len-1)/2
  pair<-rep(0,lenn)
  k=1
  for(i in 1:len){
    for(j in 1:len){
      if(i>=j){
        next
      }
      pair[k]=abs(x[i]-x[j])
      k=k+1
    }
  }
  a=median(pair)
  return(a)
}
### function for Qn
Qn<-function(x){
  len<-length(x)
  lenn<-len*(len-1)/2
  pair<-rep(0,lenn)
  k=1
  for(i in 1:len){
    for(j in 1:len){
      if(i>=j){
        next
      }
      pair[k]=abs(x[i]-x[j])
      k=k+1
    }
  }
  h=(len/2+1)*(len/2)/2
  a=sort(pair,decreasing = FALSE)
  q=a[h]
  return(q)
}
### function for LMS
LMS <- function(x)
{
  n = NROW(x);
  h = floor(n/2)+1;
  distances <- rep(0,n+1 - h);
  x <- sort(x);
  
  for(i in 1:(n+1-h))
  {
    distances[i] <- abs(x[i+h-1]-x[i]);
  }
  return(min(distances))
}

#Function which calculates average deviation
avedev <- function(x)
{
  aveX <- mean(x, na.rm = TRUE)
  return(mean(abs(x - aveX)))
}

### function for the five scale Estimators
ScaleEstimates<-function(x){
  Estimate<-rep(0.0,6)
  Estimate[1]<-MeanAD(x,center = Mean)
  Estimate[2]<-mad(x,constant = 1)
  Estimate[3]<-Sn(x,constant=1)
  Estimate[4]<-BL(x)
  Estimate[5]<-Qn(x)
  Estimate[6]<-LMS(x)
  return(Estimate)
}

#====================================#
# Sensitivity Gussian distribution   #
#====================================#

### step one: to get consistency factors
ScaleEstimates.SD<-function(x){
  Estimate<-rep(0.0,6)
  Estimate[1]<-MeanAD(x,center = Mean)
  Estimate[2]<-mad(x,constant = 1)
  Estimate[3]<-Sn(x,constant=1)
  Estimate[4]<-BL(x)
  Estimate[5]<-Qn(x)
  Estimate[6]<-LMS(x)
  Estimate[7]<-sd(x)
  return(Estimate)
}#includes the SD 

nTrial<-500
nSample<-1000
k=1
Estimates=matrix(0,0,nrow=nTrial,ncol=7)
colnames(Estimates)=c("A.Deviation","MAD","Sn","BL","Qn","LMS","SD")
for(i in 1:nTrial){
  x<-rnorm(n=nSample,mean=0,sd=1)
  Estimates[i,]<-ScaleEstimates.SD(x)
  k=k+1
  print(k)
}

boxplot(Estimates)
consFactors<-mean(Estimates[,7])/colMeans(Estimates)
consFactors_Gussian<-rep(0.0,7)
for(i in 1:7){
  consFactors_Gussian[i]=consFactors[i]
}
### Step two:Get the reference values for the estimators
set.seed(seed = 1)
x<-rnorm(100,mean=0,sd=1)
refscale<-ScaleEstimates.SD(x)*consFactors_Gussian

xAdd<-seq(from=-4,to=4,by=0.1)
nPoints<-length(xAdd)
Estimates=matrix(0.0,nrow=nPoints,ncol =7 )
rownames(Estimates)<-xAdd
colnames(Estimates)<-c("A.Deviation","MAD","Sn","BL","Qn","LMS","SD")
### step three: draw the sensitivity curve Gussian distribution)
for(i in 1:nPoints){
  ExtendedData<-c(x,xAdd[i])
  Estimates[i,]<-(ScaleEstimates.SD(ExtendedData)*consFactors_Gussian-refscale)/(1/(length(ExtendedData)))
}
Result<-melt(Estimates)
plot1<-ggplot(data=Result,aes(x=Var1,y=value,colour=as.factor(Var2)))+geom_line(size=1)+xlab("value of addded Point")+ylab("SC")+theme_bw()+
  scale_colour_manual(name="Estimator",values= c("steelblue4", "purple", "black", "lightgoldenrod4", "lightpink1", "seagreen1", " red")) +ggtitle("Gaussian model with 1 run")
plot1


#====================================#
#      Cauchy Distribution           #
#====================================#

### step one: to get consistency factors 
nSample<-1000
nTrial<-500
k=1
Estimates=matrix(0,0,nrow=nTrial,ncol=6)
colnames(Estimates)=c("A.Deviation","MAD","Sn","BL","Qn","LMS")
for(i in 1:nTrial){
  x<-rcauchy(n=nSample,location=0,scale=1)
  Estimates[i,]<-ScaleEstimates(x)
  k=k+1
  print(k)
}

boxplot(Estimates)
consFactors_Cauchy<-1/colMeans(Estimates)
consFactors_Cauchy
#consFactors_Cauchy[1]<-1 #Zach cons factor
### step two:get the reference values for the estimators
set.seed(seed = 1)
x<-rcauchy(100,location=0,scale=1)
refscale<-ScaleEstimates(x)*consFactors_Cauchy

xAdd<-seq(from=-4,to=4,by=0.1)
nPoints<-length(xAdd)
Estimates=matrix(0.0,nrow=nPoints,ncol =6 )
rownames(Estimates)<-xAdd
colnames(Estimates)<-c("A.Deviation","MAD","Sn","BL","Qn","LMS")
### step three: draw the sensitivity curve
for(i in 1:nPoints){
  ExtendedData<-c(x,xAdd[i])
  Estimates[i,]<-(ScaleEstimates(ExtendedData)*consFactors_Cauchy-refscale)/(1/length(ExtendedData))
}
Result2<-melt(Estimates)
plot2<-ggplot(data=Result2,aes(x=Var1,y=value,colour=as.factor(Var2)))+geom_line(size=1)+xlab("value of addded Point")+ylab("SC")+theme_bw()+
  scale_colour_manual(name="Estimator",values=c("steelblue4", "purple", "black", "lightgoldenrod4", "lightpink1", " red")) + ggtitle("Cauchy model with 1 run")
plot2

#==================================================#
#  Sensitiviy curve for exponential Distribuiton   #
#==================================================#

### step one: to get consistent factors 
nTrial<-500
nSample<-1000
k=1
Estimates=matrix(0,0,nrow=nTrial,ncol=7)
colnames(Estimates)=c("A.Deviation","MAD","Sn","BL","Qn","LMS","SD")
for(i in 1:nTrial){
  x<-rexp(n=nSample,r=1)
  Estimates[i,]<-ScaleEstimates.SD(x)
  k=k+1
  print(k)
}

boxplot(Estimates)
consFactors_exponential<-1/colMeans(Estimates) # scale parameter is equal to one
consFactors_exponential

### step two: Get the reference values for the estimators
set.seed(seed = 1)
x<-rexp(100,r=1)
refscale<-ScaleEstimates.SD(x)*c(consFactors_exponential)

xAdd<-seq(from=-4,to=4,by=0.1)
nPoints<-length(xAdd)
Estimates=matrix(0.0,nrow=nPoints,ncol =7 )
rownames(Estimates)<-xAdd
colnames(Estimates)<-c("A.Deviation","MAD","Sn","BL","Qn","LMS","SD")

### step three: draw the sensitivity curve 
for(i in 1:nPoints){
  ExtendedData<-c(x,xAdd[i])
  Estimates[i,]<-(ScaleEstimates.SD(ExtendedData)*consFactors_exponential-refscale)/(length(ExtendedData))
}
Result3<-melt(Estimates)
plot3<-ggplot(data=Result3,aes(x=Var1,y=value,colour=as.factor(Var2)))+geom_line(size=1)+xlab("value of addded Point")+ylab("SC")+theme_bw()+
  scale_colour_manual(name="Estimator",values=c("steelblue4", "purple", "black", "lightgoldenrod4", "lightpink1", "seagreen1", " red"))+ggtitle("Exponential model with 1 run")
plot3

#=============#
# Question 3  #
#=============#

#=====================================#
# 500 trials with Normal distribution #
#=====================================#
size=100
m=500 

range1<-seq(from=-4,to=4,by=0.1)
nPoints1<-length(range1)
MultiEstimates1 <- array(0.0, dim = c(m, nrow = nPoints1, ncol= 7))

for(j in 1:m){
  if( (j%%10)==0 ){ print(j)}
  Data1 <- rnorm(size, mean = 0, sd = 1)
  refscale1<-ScaleEstimates.SD(Data1)*consFactors_Gussian
  Estimates1 = matrix(0.0, nrow = nPoints1, ncol= 7)
  colnames(Estimates1) = c("A.Deviation","MAD","Sn","BL","Qn","LMS","SD")
  for(i in 1:nPoints1){
    ExtendedData1 <- c(Data1, range1[i])
    Estimates1[i,] <- (consFactors_Gussian*ScaleEstimates.SD(ExtendedData1) - refscale1) / (1/length(Data1))
  }
  MultiEstimates1[j,,] <- Estimates1
}

#Now calculate the mean at each location
meanEstimates <- matrix(0.0, nrow = nPoints1, ncol= 7)
colnames(meanEstimates) = c("A.Deviation","MAD","Sn","BL","Qn","LMS","SD")
rownames(meanEstimates) <- range1
for(i in 1:nPoints1){
  meanEstimates[i,] <- colMeans(MultiEstimates1[,i,])
}

#Plot the result
Result1 <- melt(meanEstimates)
plot1.1 <- ggplot(data = Result1, aes(x=Var1, y=value, colour = Var2)) +
  geom_line(size=1) +
  xlab("Value of added Point") +
  ylab("SC") +
  theme_bw() +
  ggtitle("Gaussian model with 500 runs")+
  scale_colour_manual(name = "Estimator",
                      values = c("steelblue4", "purple", "black", "lightgoldenrod4", "lightpink1", "seagreen1", " red"))
plot1.1


plots.norm <- ggarrange(plot1, plot1.1,ncol = 2, nrow = 1)
plots.norm

#=====================================#
# 500 trials with Cauchy distribution #
#=====================================#
range2<-seq(from=-4,to=4,by=0.1)
nPoints2<-length(range2)
size=100
m=500

MultiEstimates2 <- array(0.0, dim = c(m, nrow = nPoints2, ncol= 6))

for(j in 1:m){
  if( (i%%10)==0 ){ print(i)}
  Data2 <- rcauchy(size,location=0,scale=1)
  refscale2<-ScaleEstimates(Data2)*consFactors_Cauchy
  Estimates2 = matrix(0.0, nrow = nPoints2, ncol= 6)
  colnames(Estimates2) = c("A.Deviation","MAD","Sn","BL","Qn","LMS")
  for(i in 1:nPoints2){
    ExtendedData2 <- c(Data2, range2[i])
    Estimates2[i,] <- (consFactors_Cauchy*ScaleEstimates(ExtendedData2) - refscale2) / (1/length(Data2))
  }
  MultiEstimates2[j,,] <- Estimates2
}

#Now calculate the mean at each location
meanEstimates2 <- matrix(0.0, nrow = nPoints2, ncol= 6)
colnames(meanEstimates2) = c("A.Deviation","MAD","Sn","BL","Qn","LMS")
rownames(meanEstimates2) <- range2
for(i in 1:nPoints2){
  meanEstimates2[i,] <- colMeans(MultiEstimates2[,i,])
}

#Plot the result
Result2.1 <- melt(meanEstimates2)
Plot2.1 <- ggplot(data = Result2.1, aes(x=Var1, y=value, colour = Var2)) +
  geom_line(size=1) +
  xlab("Value of added Point") +
  ylab("SC") +
  theme_bw() +
  ggtitle("Cauchy model with 500 runs")+
  scale_colour_manual(name = "Estimator",
                      values = c("steelblue4", "purple", "black", "lightgoldenrod4", "lightpink1", " red"))
Plot2.1

Plots.Cauchy <- ggarrange(plot2, Plot2.1, ncol = 2, nrow = 1)
Plots.Cauchy

#===========================================#
#  500 trials with Exponential distribution # 
#===========================================#
size=100
m=500
range3<-seq(from=-4,to=4,by=0.1)
nPoints3<-length(range3)
MultiEstimates3 <- array(0.0, dim = c(m, nrow = nPoints3, ncol= 7))

for(j in 1:m){
  if( (i%%10)==0 ){ print(i)}
  Data3 <- rexp(size,r=1)
  refscale3<-ScaleEstimates.SD(Data3)*consFactors_exponential
  Estimates3 = matrix(0.0, nrow = nPoints3, ncol= 7)
  colnames(Estimates3) = c("A.Deviation","MAD","Sn","BL","Qn","LMS","SD")
  for(i in 1:nPoints3){
    ExtendedData3 <- c(Data3, range3[i])
    Estimates3[i,] <- (consFactors_exponential*ScaleEstimates.SD(ExtendedData3) - refscale3) / (1/length(Data3))
  }
  MultiEstimates3[j,,] <- Estimates3
}

#Now calculate the mean at each location
meanEstimates3 <- matrix(0.0, nrow = nPoints3, ncol= 7)
colnames(meanEstimates3) = c("A.Deviation","MAD","Sn","BL","Qn","LMS","SD")
rownames(meanEstimates3) <- range3
for(i in 1:nPoints3){
  meanEstimates3[i,] <- colMeans(MultiEstimates3[,i,])
}

#Plot the result
Result3.1 <- melt(meanEstimates3)
Plot3.1 <- ggplot(data = Result3.1, aes(x=Var1, y=value, colour = Var2)) +
  geom_line(size=1) +
  xlab("Value of added Point") +
  ylab("SC") +
  theme_bw() +
  ggtitle("Exponential model with 500 runs")+
  scale_colour_manual(name = "Estimator",
                      values = c("steelblue4", "purple", "black", "lightgoldenrod4", "lightpink1", "seagreen1", " red"))
Plot3.1

Plots.Exp <- ggarrange(plot3, Plot3.1, ncol = 2, nrow = 1)

Plots.Exp

#=================#
# Question 4      #
#=================#
rm(list = ls())

## Functions
### function for Bickel-Lehmann estimator
BL<-function(x){
  len<-length(x)
  lenn<-len*(len-1)/2
  pair<-rep(0,lenn)
  k=1
  for(i in 1:len){
    for(j in 1:len){
      if(i>=j){
        next
      }
      pair[k]=abs(x[i]-x[j])
      k=k+1
    }
  }
  a=median(pair)
  return(a)
}
### function for Qn
Qn<-function(x){
  len<-length(x)
  lenn<-len*(len-1)/2
  pair<-rep(0,lenn)
  k=1
  for(i in 1:len){
    for(j in 1:len){
      if(i>=j){
        next
      }
      pair[k]=abs(x[i]-x[j])
      k=k+1
    }
  }
  h=(len/2+1)*(len/2)/2
  a=sort(pair,decreasing = FALSE)
  q=a[h]
  return(q)
}
### function for LMS
LMS <- function(x)
{
  n = NROW(x);
  h = floor(n/2)+1;
  distances <- rep(0,n+1 - h);
  x <- sort(x);
  
  for(i in 1:(n+1-h))
  {
    distances[i] <- abs(x[i+h-1]-x[i]);
  }
  return(min(distances))
}

#Function which calculates average deviation
avedev <- function(x)
{
  aveX <- mean(x, na.rm = TRUE)
  return(mean(abs(x - aveX)))
}

#For testing purposes, normal data with two very large outliers
simdata = append(rnorm(100),rnorm(2,mean = 1000))

#To save data
pwd <- getwd()
pwd


#Parameter vector
params = list(distr = rnorm,
              avP = 1.2533,
              madP = 1.4826,
              SnP = 1.1926,
              BLP = 1.0483,
              QnP = 2.21914,
              LMSP = 0.7413)

paramsC = list(distr = rcauchy,
               avP = 1,
               madP = 1/qcauchy(0.75),
               SnP = .7071,
               BLP = 1/qcauchy(3/4, scale = 2),
               QnP = 1/qcauchy(5/8, scale = 2),
               LMSP = 1/(qcauchy(3/4, scale = 2)))  #Smallest length which contains half of the data

tableCreator <- function(n,m,params)
{
  results <- matrix(0, ncol = 6, nrow = m);
  for (i in 1:m)
  {
    data = params$distr(n)
    results[i,1] <- params$avP*avedev(data);
    results[i,2] <- params$madP*mad(data,constant = 1);
    results[i,3] <- params$SnP*Sn(data, constant = 1);
    results[i,4] <- params$BLP*BL(data);
    results[i,5] <- params$QnP*Qn(data);
    results[i,6] <- params$LMSP*LMS(data);
  }
  return(results)
}
contaminated <- function(n,epsilon,a)
{
  nClean = floor(n*(1-epsilon));
  nDirty = n - nClean;
  return(append(rnorm(nClean),rnorm(nDirty,mean = a)))
}

m <- 500;
nlist <- c(10,20,40,60,80,100,200,1000);
resultsAV <- matrix(0,nrow = 8, ncol = 6);
resultsSD <- matrix(0,nrow = 8, ncol = 6);
p <- 1;
params <- paramsC;
for (n in nlist)
{
  print(n)
  M = tableCreator(n,m,params)
  
  if(identical(rnorm,params$distr))
  {
    print("saving gaussian data")
    filename <- paste("n",n,"m",m,"normal.txt")
  } else if (identical(rcauchy, params$distr))
  {
    print("saving cauchy data")
    filename <- paste("n",n,"m",m,"cauchy.txt")
  }
  
  write.table(M,filename,sep = ' ')
  resultsAV[p,] <- apply(M,2,mean)
  resultsSD[p,] <- apply(M,2,var)*n/(resultsAV[p,]^2)
  p <- p + 1;
}

resultsAV.df <- data.frame(cbind(nlist,resultsAV))
colnames(resultsAV.df) <- c("N", "AvDev", "MAD", "Sn", "Bickel-Lehman","Qn","LMS")

resultsSD.df <- data.frame(cbind(nlist, resultsSD))
colnames(resultsSD.df) <- c("N", "AvDev", "MAD", "Sn", "Bickel-Lehman","Qn","LMS")

if(identical(rnorm,params$distr))
{
  print("saving gaussian data")
  filename1 <- paste("resultsAV_n",n,"_m_",m,"_normal.Rda")
  filename2 <- paste("resultsSD_n",n,"_m_",m,"_normal.Rda")
} else if (identical(rcauchy, params$distr))
{
  print("saving cauchy data")
  filename1 <- paste("resultsAV_n",n,"_m_",m,"_cauchy.Rda")
  filename2 <- paste("resultsSD_n",n,"_m_",m,"_cauchy.Rda")
}

save(resultsAV.df,file = filename1)
save(resultsSD.df, file = filename2)

#===============#
# Question 5    #
#===============#

n = 100;
a = 10;
m = 500;

params = list(distr = rnorm,
              avP = 1.2533,
              madP = 1.4826,
              SnP = 1.1926,
              BLP = 1.0483,
              QnP = 2.21914,
              LMSP = 0.7413)
###
epsilonList <- seq(from = 0.01, to = 0.49, by = 0.01);
plot1 <- matrix(0, nrow = length(epsilonList), ncol = 6);
p <- 1;

for(epsilon in epsilonList)
{
  params$distr <- (function(x){contaminated(x,epsilon,a)});
  M <- tableCreator(n,m,params)
  plot1[p,] <- apply(M, 2, mean);
  p <- p+1;
}
contaminated
epsilonPlot <- data.frame(plot1)
colnames(epsilonPlot) <- c("AverageDeviation", "MAD", "Sn", "BL", "Qn","LMS")
epsilonPlot['Epsilon'] <- epsilonList

### 

epsilon = 0.2
aList <- seq(from = 0, to = 10, by = 0.1)
plot2 <- matrix(0,nrow = length(aList), ncol = 6)
p <- 1;

for(a in aList)
{
  params$distr <- (function(x){contaminated(x,epsilon,a)});
  M <- tableCreator(n,m,params)
  plot2[p,] <- apply(M, 2, mean);
  p <- p+1;
}

aPlot <- data.frame(plot2)
colnames(aPlot) <- c("AverageDeviation", "MAD", "Sn", "BL", "Qn","LMS")
aPlot['a'] <- aList;

library(ggplot2)

ggplot(epsilonPlot, aes(x = Epsilon)) + geom_line(aes(y = AverageDeviation, color = 'Av. Dev.')) + geom_line(aes(y = MAD, color = 'MADN')) + 
  geom_line(aes(y = Sn, color = 'Sn')) + geom_line(aes(y = BL, color = 'BL')) + geom_line(aes(y = Qn, color = 'Qn')) + 
  geom_line(aes(y = LMS, color = 'LMS')) + ylab('Mean Scale Estimator') + ggtitle("Effect of contamination on each estimate") + 
  scale_color_manual(name = "Legend", values = c('Av. Dev.' = 'red', MADN = 'orange', Sn = 'green', BL = 'blue', Qn = 'purple', LMS = 'violet'))


ggplot(aPlot, aes(x = a)) + geom_line(aes(y = AverageDeviation, color = 'Av. Dev.')) + geom_line(aes(y = MAD, color = 'MADN')) + 
  geom_line(aes(y = Sn, color = 'Sn')) + geom_line(aes(y = BL, color = 'BL')) + geom_line(aes(y = Qn, color = 'Qn')) + 
  geom_line(aes(y = LMS, color = 'LMS')) + ylab('Mean Scale Estimator') + ggtitle("Effect of large deviations") + 
  scale_color_manual(name = "Legend", values = c('Av. Dev.' = 'red', MADN = 'orange', Sn = 'green', BL = 'blue', Qn = 'purple', LMS = 'violet'))


#=============#
# Question 6  #
#=============#

n = 100;
a = 10;
m = 500;

params = list(distr = rnorm,
              avP = 1.2533,
              madP = 1.4826,
              SnP = 1.1926,
              BLP = 1.0483,
              QnP = 2.21914,
              LMSP = 0.7413)
###
epsilonList1 <- seq(from = 0.01, to = 0.49, by = 0.01);
plot3 <- matrix(0, nrow = length(epsilonList1), ncol = 6);
p <- 1;

#REDEFINE CONTAMINATED
epsilon = 0.2
aa = c(-a,a)
contaminated1 <- function(n,epsilon,aa)
{
  nClean <- (n*(1-epsilon));
  nDirty <- n - nClean;
  return(append(rnorm(nClean),rnorm(nDirty,mean = aa, sd=1)));
}

### Defining plot 1
for(epsilon in epsilonList1)
{
    params$distr <- (function(x){contaminated1(x,epsilon,aa)});
    M <- tableCreator(n,m,params)
    plot3[p,] <- apply(M, 2, mean);
    p <- p+1;
}

epsilonPlot2 <- data.frame(plot3)
colnames(epsilonPlot2) <- c("AverageDeviation", "MAD", "Sn", "BL", "Qn","LMS")
epsilonPlot2['Epsilon'] <- epsilonList1

### Defining plot 2

aList1 <- seq(from = 0, to = 10, by = 0.1)
plot4 <- matrix(0,nrow = length(aList1), ncol = 6)
p <- 1;

for(a in aList1)
{
  aa=c(-a,a)
  params$distr <- (function(x){contaminated1(x,epsilon,aa)});
  M <- tableCreator(n,m,params)
  plot4[p,] <- apply(M, 2, mean);
  p <- p+1;
}

aPlot2 <- data.frame(plot4)
colnames(aPlot2) <- c("AverageDeviation", "MAD", "Sn", "BL", "Qn","LMS")
aPlot2['a'] <- aList1;

### plot 
library(ggplot2)

ggplot(epsilonPlot2, aes(x = Epsilon)) + geom_line(aes(y = AverageDeviation, color = 'Av. Dev.')) + geom_line(aes(y = MAD, color = 'MADN')) + 
  geom_line(aes(y = Sn, color = 'Sn')) + geom_line(aes(y = BL, color = 'BL')) + geom_line(aes(y = Qn, color = 'Qn')) + 
  geom_line(aes(y = LMS, color = 'LMS')) + ylab('Mean Scale Estimator') + ggtitle("Effect of contamination on each estimate") + 
  scale_color_manual(name = "Legend", values = c('Av. Dev.' = 'red', MADN = 'orange', Sn = 'green', BL = 'blue', Qn = 'purple', LMS = 'violet'))


ggplot(aPlot2, aes(x = a)) + geom_line(aes(y = AverageDeviation, color = 'Av. Dev.')) + geom_line(aes(y = MAD, color = 'MADN')) + 
  geom_line(aes(y = Sn, color = 'Sn')) + geom_line(aes(y = BL, color = 'BL')) + geom_line(aes(y = Qn, color = 'Qn')) + 
  geom_line(aes(y = LMS, color = 'LMS')) + ylab('Mean Scale Estimator') + ggtitle("Effect of large deviations") + 
  scale_color_manual(name = "Legend", values = c('Av. Dev.' = 'red', MADN = 'orange', Sn = 'green', BL = 'blue', Qn = 'purple', LMS = 'violet'))

#==============#
# Question 7   #
#==============#
diabetes <- read.csv("H:/Masters_of_Statistics/Semester_4/Robust Statistics/Project/Data set/diabetes.csv")
head(diabetes)

#Takes M trials and a sample of size n from dataset X with given parameters
scale.estimates <- function(data,params)
{
  results <- rep(0,7)
  results[1] <- params$avP*MeanAD(data)#avedev(data);
  results[2] <- params$madP*mad(data,constant = 1);
  results[3] <- params$SnP*robustbase::Sn(data, constant = 1);
  results[4] <- params$BLP*BL(data);
  results[5] <- params$QnP*robustbase::Qn(data, constant = 1);
  results[6] <- params$LMSP*LMS(data);
  results[7] <- sd(data);
  return(results)
}

insulin <- diabetes$Insulin[diabetes$Insulin != 0]
bmi <- diabetes$BMI[diabetes$BMI != 0]

data <- bmi

estNames <-  c("meanAD", "MAD","Sn","BL","Qn","LMS","SD")
scale <- scale.estimates(data, params)
p <- 1
outliers.df <- (data - median(data))/scale[1]
for(i in 2:length(scale))
{
  outliers.df <- cbind(outliers.df, (data - median(data))/scale[i])
}

outliers.df <- abs(outliers.df)
colnames(outliers.df) <- estNames
outliers.df <- data.frame(outliers.df)

diabet.class <- diabetes%>%
  filter(BMI != 0)%>%
  select(BMI, Outcome)

mean(diabet.class[outliers.df$meanAD > 3,]$Outcome)
mean(diabet.class[outliers.df$MAD > 3, ]$Outcome)
mean(diabet.class[outliers.df$Sn > 3, ]$Outcome)
mean(diabet.class[outliers.df$BL > 3, ]$Outcome)
mean(diabet.class[outliers.df$Qn > 3, ]$Outcome)
mean(diabet.class[outliers.df$LMS > 3, ]$Outcome)
mean(diabet.class[outliers.df$SD > 3, ]$Outcome)

library(ggplot2)
outliers.df$rowID = 1:nrow(outliers.df)
ggplot(outliers.df, aes(x = rowID, y = meanAD, color = as.factor(meanAD > 3))) + 
  geom_point() + 
  geom_hline(yintercept =  3, color = 'red') + 
  scale_color_manual(values = c("black", 'red'))

ggplot(outliers.df, aes(x = rowID, y = MAD, color = as.factor(MAD > 3))) + 
  geom_point() + 
  geom_hline(yintercept =  3, color = 'red') + 
  scale_color_manual(values = c("black", 'red'))

ggplot(outliers.df, aes(x = rowID, y = Sn, color = as.factor(Sn > 3))) + 
  geom_point() + 
  geom_hline(yintercept =  3, color = 'red') + xlab("Index")+
  scale_color_manual(values = c("black", 'red'))

ggplot(outliers.df, aes(x = rowID, y = Qn, color = as.factor(Qn > 3))) + 
  geom_point() + 
  geom_hline(yintercept =  3, color = 'red') + 
  scale_color_manual(values = c("black", 'red'))

ggplot(outliers.df, aes(x = rowID, y = BL, color = as.factor(BL > 3))) + 
  geom_point() + 
  geom_hline(yintercept =  3, color = 'red') + 
  scale_color_manual(values = c("black", 'red'))

ggplot(outliers.df, aes(x = rowID, y = LMS, color = as.factor(LMS > 3))) + 
  geom_point() + 
  geom_hline(yintercept =  3, color = 'red') + 
  scale_color_manual(values = c("black", 'red'))

ggplot(outliers.df, aes(x = rowID, y = SD, color = as.factor(SD > 3))) + 
  geom_point() + 
  geom_hline(yintercept =  3, color = 'red') + 
  scale_color_manual(values = c("black", 'red'))


results <- rbind( avdev = scale.estimates(data[outliers.df$meanAD <= 3], params),
                  MAD = scale.estimates(data[outliers.df$MAD <= 3], params),
                  Sn = scale.estimates(data[outliers.df$Sn <= 3], params),
                  BL = scale.estimates(data[outliers.df$BL <= 3], params),
                  Qn = scale.estimates(data[outliers.df$Qn <= 3], params),
                  LMS = scale.estimates(data[outliers.df$LMS <= 3], params),
                  SD = scale.estimates(data[outliers.df$SD <= 3], params),
                  None = scale)


colnames(results) <- estNames
results

