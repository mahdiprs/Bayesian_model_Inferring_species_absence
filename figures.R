#Creates plots found in paper
#Title: Inferring species absence from zero-sighting records using
#       analytical Bayesian models with population growth
#Authors: B.Barnes, F.Giannini, M.Parsa and D.Ramsey
 
library(reshape2)
library(ggplot2)
library(RColorBrewer)

#FUNCTIONS
readAnderson<-function(fileNameEx)
{
  #Read in Anderson model ouput data 
  npvToErad_SII<-read.csv(paste("Data\\npvToErad_SII_", fileNameEx, ".csv", sep=""), header=T)
  npvToErad_SII<-as.matrix(npvToErad_SII[,-1])
  yrToErad_SII<-read.csv(paste("Data\\yrToErad_SII_", fileNameEx, ".csv", sep=""), header=T)
  yrToErad_SII<-as.matrix(yrToErad_SII[,-1])
  pofOptim_SII<-read.csv(paste("Data\\pofOptim_SII_", fileNameEx, ".csv", sep=""), header=T)
  pofOptim_SII<-as.matrix(pofOptim_SII[,-1])
  annSurveyCost_SII<-read.csv(paste("Data\\annSurveyCost_SII_", fileNameEx, ".csv", sep=""), header=T)
  annSurveyCost_SII<-as.matrix(annSurveyCost_SII[,-1])
  
  melted_cost <- melt(t(npvToErad_SII))
  melted_years <- melt(t(yrToErad_SII))
  melted_pof <-melt(t(pofOptim_SII))
  melted_surveyCost<-melt(t(annSurveyCost_SII))
  
  allData<-cbind(melted_cost, melted_years$value, melted_pof$value, melted_surveyCost$value)
  allData$Var1<-pdSURange_SII[allData$Var1]
  allData$Var2<-perCovRange_SII[allData$Var2]
  
  names(allData)<-c("pd", "prp", "TotalCost", "nSurveys", "pof", "annSurveyCost")
  completeData<- allData[complete.cases(allData), ]
  return(completeData)
} 

numInfestedUnits<-function(Z0, r, K, kappaV)
{
  expectedZ<-rep(0, length(kappaV)+1)
  expectedZ[1]<-Z0
  
  for(t in 1:length(kappaV))
  {  
   expectedZ[t+1]<-expectedZ[t]+r*expectedZ[t]*(K-expectedZ[t])/K
  }
  
  return(expectedZ)  
}

numSUSearched<-function(probP, pofT, Z, numSurveys, N, gammaV)
{
  H<-1/sum(Z[2:(numSurveys+1)])
  nSU<-(1-(((1-probP)/probP)*((1-pofT)/pofT))^H)*(N/gammaV)
  return(nSU)
} 

numSUSearched_rApprox<-function(probP, pofT, r, numSurveys, N, gammaV)
{
  H_approx<-r/((1+r)^(numSurveys+1)-(1+r))
  nSU<-(1-(((1-probP)/probP)*((1-pofT)/pofT))^H_approx)*(N/gammaV)
  return(nSU)
} 

cost<-function(nSurveys, n, N, gammaV, errorCost)
{
  totalCost<-nSurveys*(n/N)*(N*log(1-gammaV)/log(1-0.01))+errorCost
  return(totalCost)
}
  
solveForDelta<-function(lambda, p, T, kappa)
{  
  delta<-seq(0.0000001, 0.999999, by=0.000001)
  s=0
  solution<-rep(NA, length(kappa))
  
  G<-exp(lambda*((1-delta)*s-1))
  H<-exp(lambda*((1-delta)-1))
  
  f1<-(1-p+p*G)/(1-p+p*H)
  
  solution[1]<-delta[which(f1>T)[1]]
  
  for(i in 2:length(kappa))
  {
    G<-exp(lambda*((1-delta)*G-1))
    H<-exp(lambda*((1-delta)*H-1))
    fn<-(1-p+p*G)/(1-p+p*H)
    
    solution[i]<-delta[which(fn>T)[1]]
  }
  
  return(cbind(kappa, solution))
}



#PARAMETERS
#Barnes equation parameters 
pof=0.95 #prob. of freedom threshold
z0<-1 #initial no. of infested units
p<-0.75 #prob. of presence
r<-1 #deterministic growth 
K<-1000 #carrying capacity 
lambda<-2

#Anderson app parameters
N<-5000 #no. of SUs in a MZ
nV<-500 #no. of values in range from Anderson model
Cm_star<-(400+1)*209897 #Re-control cost from Anderson
errorCost<-Cm_star*(1-0.9561)

gammaVector<-seq(0.01,0.99,length.out = nV)
kappaVector<-1:15 #no. of surveys

pdSURange_SII<-seq(0.01, 0.99, length.out = nV)
perCovRange_SII<-seq(0.01, 0.99, length.out = nV)

#Plot parameters
plotSurveys_I<-c(3, 7, 11, 15) #for no growth, det. growth plots
plotSurveys_II<-c(3, 5, 7, 9) #for det. vs stochastic growth plots
legText_I<-expression(paste( kappa, "=3 ", ".... ", kappa, "=7 ",".-.- ", 
                           kappa, "=11 ", "- - ", kappa, "=15"))
legText_II<-expression(paste( kappa, "=3 ", ".... ", kappa, "=5 ",".-.- ", 
                           kappa, "=7 ", "- - ", kappa, "=9"))
colA<-"skyblue"
colB<-"black"
cx<- 1.4 # legend txt
n<-matrix(NA, nrow=length(plotSurveys_I), ncol=nV)
n_det<-matrix(NA, nrow=length(plotSurveys_II), ncol=nV)
n_stoch<-matrix(NA, nrow=length(plotSurveys_II), ncol=nV)

#PLOT
#Fig 1 (a) and (b)
#No growth, no spread, 1 MZ

#Read in Anderson model ouput data
andersonData<-readAnderson("noGrowth")

extractCurve<-andersonData[which(andersonData$nSurveys==plotSurveys_I[1]),]
extractCurve<-extractCurve[which(extractCurve$pof<0.951),]

setEPS()
postscript("Fig1a_AreaSearched_NoGrowth.eps") 
par(mai=c(2.5,0.8,0.4,0.8)+0.1, xpd=TRUE)

plot(extractCurve$pd, extractCurve$prp, type='l', lty=1, lwd=4, col=colA, 
     xlim=c(0,1), ylim=c(0,1), xlab=expression(paste("search-effort (", gamma, " )", sep='')), 
     ylab="proportion of area searched (n/N)",
     xaxs='i', yaxs='i',  xaxt="n", yaxt="n",cex.lab=cx)
labelXY<-seq(0,1,0.2)
labelXY[1]<-0
labelXY[length(labelXY)]<-1
axis(side=1, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)
axis(side=2, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)

for(i in 2:length(plotSurveys_I))
{  
  extractCurve<-andersonData[which(andersonData$nSurveys==plotSurveys_I[i]),]
  extractCurve<-extractCurve[which(extractCurve$pof<0.951),]
  lines(extractCurve$pd, extractCurve$prp, lty=1, col=colA, lwd=4)
} 

z<-numInfestedUnits(z0, 0, 1, kappaVector) #no growth

for(i in 1:length(plotSurveys_I))
{  
  n[i,]<-numSUSearched(p, pof, z, plotSurveys_I[i], N, gammaVector)
  barnesCurve<-cbind(gammaVector, n[i,])
  if(length(which(barnesCurve[,2]>N))!=0)
  {  
    barnesCurve<-barnesCurve[-which(barnesCurve[,2]>N),] #remove rows where n >N
  } 
  
  barnesCurve[,2]<-barnesCurve[,2]/N #change to prp
  lines(barnesCurve[,1], barnesCurve[,2], lty=i+1, lwd=2, col=colB) #lty=i
}


legend("topright", legend=legText_I, lty=2, lwd=1, bty='n', inset=c(0, -0.15),cex=cx) 
dev.off()


setEPS()
postscript("Fig1b_Cost_NoGrowth.eps")
par(mai=c(2.5,0.8,0.4,0.8)+0.1, xpd=TRUE)

extractCurve<-andersonData[which(andersonData$nSurveys==plotSurveys_I[1]),]
extractCurve<-extractCurve[which(extractCurve$pof<0.951),]

plot(extractCurve$pd, extractCurve$TotalCost/10^6, type='l', lty=1, lwd=4, col=colA, 
     xlim=c(0,1), ylim=c(5,12), xlab=expression(paste("search-effort (", gamma, " )", sep='')), 
     ylab=expression('total expected cost (x10'^"6"*')'),
     xaxs='i', yaxs='i',  xaxt="n", yaxt="n",cex.lab=cx)
labelXY<-seq(0,1,0.2)
labelXY[1]<-0
labelXY[length(labelXY)]<-1
axis(side=1, at=seq(0,1,0.2), labels=labelXY,cex.axis=1.2)
axis(side=2, at=seq(5, 12, 1),cex.axis=1.2)

for(i in 2:length(plotSurveys_I))
{  
  extractCurve<-andersonData[which(andersonData$nSurveys==plotSurveys_I[i]),]
  extractCurve<-extractCurve[which(extractCurve$pof<0.951),]
  lines(extractCurve$pd, extractCurve$TotalCost/10^6, lty=1, col=colA, lwd=4)
} 

for(i in 1:length(plotSurveys_I))
{  
  totalCost<-cost(plotSurveys_I[i], n[i,], N, gammaVector, errorCost)
  barnesCurve<-cbind(gammaVector, n[i,], totalCost)
  
  if(length(which(barnesCurve[,2]>N))!=0)
  {  
    barnesCurve<-barnesCurve[-which(barnesCurve[,2]>N),] #remove rows where n >N
  }
  
  lines(barnesCurve[,1], barnesCurve[,3]/10^6, lty=i+1, lwd=2, col=colB)
}

legend("topright", legend=legText_I, lty=2, lwd=2, bty='n', inset=c(0, -0.15),cex=cx) 
dev.off()


#Fig 1 (c) and (d)
#Growth, no spread, 1 MZ

#Read in Anderson model ouput data 
andersonData<-readAnderson("Growth")

extractCurve<-andersonData[which(andersonData$nSurveys==plotSurveys_I[1]),]
extractCurve<-extractCurve[which(extractCurve$pof<0.951),]

setEPS()
postscript("Fig1c_AreaSearched_Growth.eps")
par(mai=c(2.5,0.8,0.4,0.8)+0.1, xpd=TRUE)

plot(extractCurve$pd, extractCurve$prp, type='l', lty=1, lwd=4, col=colA, 
     xlim=c(0,1), ylim=c(0,1), xlab=expression(paste("search-effort (", gamma, " )", sep='')), 
     ylab="proportion of area searched (n/N)",
     xaxs='i', yaxs='i',  xaxt="n", yaxt="n", cex.lab=cx)
labelXY<-seq(0,1,0.2)
labelXY[1]<-0
labelXY[length(labelXY)]<-1
axis(side=1, at=seq(0,1,0.2), labels=labelXY,cex.axis=1.2)
axis(side=2, at=seq(0,1,0.2), labels=labelXY,cex.axis=1.2)

for(i in 2:length(plotSurveys_I))
{  
  extractCurve<-andersonData[which(andersonData$nSurveys==plotSurveys_I[i]),]
  extractCurve<-extractCurve[which(extractCurve$pof<0.951),]
  lines(extractCurve$pd, extractCurve$prp, lty=1, lwd=4, col=colA)
} 

z<-numInfestedUnits(z0, r, K, kappaVector)

for(i in 1:length(plotSurveys_I))
{
  n[i,]<-numSUSearched(p, pof, z, plotSurveys_I[i], N, gammaVector)
  barnesCurve<-cbind(gammaVector, n[i,])

  if(length(which(barnesCurve[,2]>N)!=0))
  {
   barnesCurve<-barnesCurve[-which(barnesCurve[,2]>N),] #remove rows where n >N
  }
    
  barnesCurve[,2]<-barnesCurve[,2]/N #change to prp
  lines(barnesCurve[,1], barnesCurve[,2], lty=i+1, lwd=2, col=colB)
}

legend("topright", legend=legText_I, lty=2, lwd=1, bty='n',inset=c(0, -0.15),cex=cx)
dev.off()

extractCurve<-andersonData[which(andersonData$nSurveys==plotSurveys_I[1]),]
extractCurve<-extractCurve[which(extractCurve$pof<0.951),]

setEPS()
postscript("Fig1d_Cost_Growth.eps")
par(mai=c(2.5,0.8,0.4,0.8)+0.1, xpd=TRUE)

plot(extractCurve$pd, extractCurve$TotalCost/10^6, type='l', lty=1, lwd=4, col=colA, 
     xlim=c(0,1), ylim=c(3.6,5.6), xlab=expression(paste("search-effort (", gamma, ")", sep='')), 
     ylab=expression('total expected cost (x10'^"6"*')'),
     xaxs='i', yaxs='i',  xaxt="n", yaxt="n", cex.lab=1.4)
labelXY<-seq(0,1,0.2)
labelXY[1]<-0
labelXY[length(labelXY)]<-1
axis(side=1, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)
axis(side=2, at=seq(3.6, 5.6, 0.2), cex.axis=1.2)

for(i in 2:length(plotSurveys_I))
{  
  extractCurve<-andersonData[which(andersonData$nSurveys==plotSurveys_I[i]),]
  extractCurve<-extractCurve[which(extractCurve$pof<0.951),]
  lines(extractCurve$pd, extractCurve$TotalCost/10^6, lty=1, col=colA, lwd=4)
} 

for(i in 1:length(plotSurveys_I))
{  
  totalCost<-cost(plotSurveys_I[i], n[i,], N, gammaVector, errorCost)
  barnesCurve<-cbind(gammaVector, n[i,], totalCost)
  
  if(length(which(barnesCurve[,2]>N))!=0)
  {  
   barnesCurve<-barnesCurve[-which(barnesCurve[,2]>N),] #remove rows where n >N
  } 
  
  lines(barnesCurve[,1], barnesCurve[,3]/10^6, lty=i+1, lwd=2, col=colB)
}

legend("topright", legend=legText_I, lty=2, lwd=2, bty='n', inset=c(0, -0.15),cex=cx) 
dev.off()

#Fig 3 (a) and (b)
#Comparison between deterministic and stochastic growth
postscript("Fig3a_AreaSearched_Stoch.eps")
par(mai=c(2.5,0.9,0.4,0.8)+0.1, xpd=TRUE)

plot(gammaVector, gammaVector, type='n', lty=2, lwd=2, col=colA, 
     xlim=c(0,1), ylim=c(0,1), xlab=expression(paste("search-effort (", gamma, ")", sep='')), 
     ylab="proportion of area searched (n/N)",
     xaxs='i', yaxs='i',  xaxt="n", yaxt="n", cex.lab=1.4)
labelXY<-seq(0,1,0.2)
labelXY[1]<-0
labelXY[length(labelXY)]<-1
axis(side=1, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)
axis(side=2, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)

for(i in 1:length(plotSurveys_II))
{  
  n_det[i,]<-numSUSearched_rApprox(p, pof, r, plotSurveys_II[i], N, gammaVector)
  detCurve<-cbind(gammaVector, n_det[i,])
  
  if(length(which(detCurve[,2]>N)!=0))
  {
    detCurve<-detCurve[-which(detCurve[,2]>N),] #remove rows where n >N
  }
  
  detCurve[,2]<-detCurve[,2]/N #change to prp
  lines(detCurve[,1], detCurve[,2], lty=i+1, lwd=2, col=colA)
  
}

stochData<-solveForDelta(lambda, p, pof, kappaVector)
stochData<-as.data.frame(stochData)
colnames(stochData)<-c("kappa", "gammaN")

for(i in 1:length(plotSurveys_II))
{  
  n_stoch[i,]<-N*stochData[which(stochData$kappa==plotSurveys_II[i]),]$gammaN/gammaVector
  stochCurve<-cbind(gammaVector, n_stoch[i,])
  
  if(length(which(stochCurve[,2]>N)!=0))
  {
    stochCurve<-stochCurve[-which(stochCurve[,2]>N),] #remove rows where n >N
  }
  
  stochCurve[,2]<-stochCurve[,2]/N #change to prp
  lines(stochCurve[,1], stochCurve[,2], lty=i+1, lwd=2, col=colB)
  
}

legend("topright", legend=legText_II, lty=2, lwd=1, bty='n',inset=c(0, -0.15), cex=cx) 

dev.off()


setEPS()
postscript("Fig3b_Cost_Stoch.eps")
par(mai=c(2.5,0.8,0.4,0.8)+0.1, xpd=TRUE)

for(i in 1:length(plotSurveys_II))
{  
  totalCost<-cost(plotSurveys_II[i], n_det[i,], N, gammaVector, errorCost)
  detCurve<-cbind(gammaVector, n_det[i,], totalCost)
  
  if(length(which(detCurve[,2]>N))!=0)
  {  
    detCurve<-detCurve[-which(detCurve[,2]>N),] #remove rows where n >N
  }
  
  if(i==1)
  {
    plot(detCurve[,1], detCurve[,3]/10^6, type='l', lty=i+1, lwd=2, col=colA, 
         xlim=c(0,1), ylim=c(3.6,6), xlab=expression(paste("search-effort (", gamma, ")", sep='')), 
         ylab=expression('total expected cost (x10'^"6"*')'),
         xaxs='i', yaxs='i',  xaxt="n", yaxt="n", cex.lab=1.4)
    labelXY<-seq(0,1,0.2)
    labelXY[1]<-0
    labelXY[length(labelXY)]<-1
    axis(side=1, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)
    axis(side=2, at=seq(3.6, 6, 0.4), cex.axis=1.2)
    
  }else{  
    lines(detCurve[,1], detCurve[,3]/10^6, lty=i+1, lwd=2, col=colA)
  }
  
}

for(i in 1:length(plotSurveys_II))
{  
  totalCost<-cost(plotSurveys_II[i], n_stoch[i,], N, gammaVector, errorCost)
  stochCurve<-cbind(gammaVector, n_stoch[i,], totalCost)
  
  if(length(which(stochCurve[,2]>N))!=0)
  {  
    stochCurve<-stochCurve[-which(stochCurve[,2]>N),] #remove rows where n >N
  }
  
  lines(stochCurve[,1], stochCurve[,3]/10^6, lty=i+1, lwd=2, col=colB)
}

legend("topright", legend=legText_II, lty=2, lwd=1, bty='n', inset=c(0, -0.15),cex=cx) 
dev.off()

#Fig 3 (c) and (d)
#Stochastic vs deterministic growth
#Different pof thresholds

postscript("Fig3c_T0.85_Stoch.eps")
par(mai=c(2.5,0.9,0.4,0.8)+0.1, xpd=TRUE)

plot(gammaVector, gammaVector, type='n', lty=2, lwd=2, col=colA, 
     xlim=c(0,1), ylim=c(0,1), xlab=expression(paste("search-effort (", gamma, ")", sep='')), 
     ylab="proportion of area searched (n/N)",
     xaxs='i', yaxs='i',  xaxt="n", yaxt="n", cex.lab=1.4)
labelXY<-seq(0,1,0.2)
labelXY[1]<-0
labelXY[length(labelXY)]<-1
axis(side=1, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)
axis(side=2, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)

for(i in 1:length(plotSurveys_II))
{  
  n_det[i,]<-numSUSearched_rApprox(p, 0.85, r, plotSurveys_II[i], N, gammaVector) #pof=0.85
  detCurve<-cbind(gammaVector, n_det[i,])
  
  if(length(which(detCurve[,2]>N)!=0))
  {
    detCurve<-detCurve[-which(detCurve[,2]>N),] #remove rows where n >N
  }
  
  detCurve[,2]<-detCurve[,2]/N #change to prp
  if(i==1){
    lines(detCurve[,1], detCurve[,2], lty=i+1, lwd=4, col=colA) 
  } else {
    lines(detCurve[,1], detCurve[,2], lty=i+1, lwd=2, col=colA)
  }  
}

stochData<-solveForDelta(lambda, p, 0.85, kappaVector) #pof=0.85
stochData<-as.data.frame(stochData)
colnames(stochData)<-c("kappa", "gammaN")

for(i in 1:length(plotSurveys_II))
{  
  n_stoch[i,]<-N*stochData[which(stochData$kappa==plotSurveys_II[i]),]$gammaN/gammaVector
  stochCurve<-cbind(gammaVector, n_stoch[i,])
  
  if(length(which(stochCurve[,2]>N)!=0))
  {
    stochCurve<-stochCurve[-which(stochCurve[,2]>N),] #remove rows where n >N
  }
  
  stochCurve[,2]<-stochCurve[,2]/N #change to prp
  lines(stochCurve[,1], stochCurve[,2], lty=i+1, lwd=2, col=colB)
  
}

legend("topright", legend=legText_II, lty=2, lwd=1, bty='n',inset=c(0, -0.15), cex=cx) 

dev.off()

postscript("Fig3d_T0.75_Stoch.eps")
par(mai=c(2.5,0.9,0.4,0.8)+0.1, xpd=TRUE)

plot(gammaVector, gammaVector, type='n', lty=2, lwd=2, col=colA, 
     xlim=c(0,1), ylim=c(0,1), xlab=expression(paste("search-effort (", gamma, ")", sep='')), 
     ylab="proportion of area searched (n/N)",
     xaxs='i', yaxs='i',  xaxt="n", yaxt="n", cex.lab=1.4)
labelXY<-seq(0,1,0.2)
labelXY[1]<-0
labelXY[length(labelXY)]<-1
axis(side=1, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)
axis(side=2, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)

for(i in 1:length(plotSurveys_II))
{  
  n_det[i,]<-numSUSearched_rApprox(p, 0.75, r, plotSurveys_II[i], N, gammaVector) #pof=0.75
  detCurve<-cbind(gammaVector, n_det[i,])
  
  if(length(which(detCurve[,2]>N)!=0))
  {
    detCurve<-detCurve[-which(detCurve[,2]>N),] #remove rows where n >N
  }
  
  detCurve[,2]<-detCurve[,2]/N #change to prp
  lines(detCurve[,1], detCurve[,2], lty=i+1, lwd=2, col=colA)
}

stochData<-solveForDelta(lambda, p, 0.75, kappaVector) #pof=0.75
stochData<-as.data.frame(stochData)
colnames(stochData)<-c("kappa", "gammaN")

for(i in 1:length(plotSurveys_II))
{  
  n_stoch[i,]<-N*stochData[which(stochData$kappa==plotSurveys_II[i]),]$gammaN/gammaVector
  stochCurve<-cbind(gammaVector, n_stoch[i,])
  
  if(length(which(stochCurve[,2]>N)!=0))
  {
    stochCurve<-stochCurve[-which(stochCurve[,2]>N),] #remove rows where n >N
  }
  
  stochCurve[,2]<-stochCurve[,2]/N #change to prp
  lines(stochCurve[,1], stochCurve[,2], lty=i+1, lwd=2, col=colB)
  
}

legend("topright", legend=legText_II, lty=2, lwd=1, bty='n',inset=c(0, -0.15), cex=cx) 

dev.off()

#Fig 3 (e) and (f)
#Stochastic vs deterministic growth
#Different growth rates
postscript("Fig3e_r0.5lambda1.5_Stoch.eps")
par(mai=c(2.5,0.9,0.4,0.8)+0.1, xpd=TRUE)

plot(gammaVector, gammaVector, type='n', lty=2, lwd=2, col=colA, 
     xlim=c(0,1), ylim=c(0,1), xlab=expression(paste("search-effort (", gamma, ")", sep='')), 
     ylab="proportion of area searched (n/N)",
     xaxs='i', yaxs='i',  xaxt="n", yaxt="n", cex.lab=1.4)
labelXY<-seq(0,1,0.2)
labelXY[1]<-0
labelXY[length(labelXY)]<-1
axis(side=1, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)
axis(side=2, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)

for(i in 1:length(plotSurveys_II))
{  
  n_det[i,]<-numSUSearched_rApprox(p, pof, 0.5, plotSurveys_II[i], N, gammaVector) #r=0.5
  detCurve<-cbind(gammaVector, n_det[i,])
  
  if(length(which(detCurve[,2]>N)!=0))
  {
    detCurve<-detCurve[-which(detCurve[,2]>N),] #remove rows where n >N
  }
  
  detCurve[,2]<-detCurve[,2]/N #change to prp
  lines(detCurve[,1], detCurve[,2], lty=i+1, lwd=2, col=colA)
}

stochData<-solveForDelta(1.5, p, pof, kappaVector) #lambda=1.5
stochData<-as.data.frame(stochData)
colnames(stochData)<-c("kappa", "gammaN")

for(i in 1:length(plotSurveys_II))
{  
  n_stoch[i,]<-N*stochData[which(stochData$kappa==plotSurveys_II[i]),]$gammaN/gammaVector
  stochCurve<-cbind(gammaVector, n_stoch[i,])
  
  if(length(which(stochCurve[,2]>N)!=0))
  {
    stochCurve<-stochCurve[-which(stochCurve[,2]>N),] #remove rows where n >N
  }
  
  stochCurve[,2]<-stochCurve[,2]/N #change to prp
  lines(stochCurve[,1], stochCurve[,2], lty=i+1, lwd=2, col=colB)
  
}

legend("topright", legend=legText_II, lty=2, lwd=1, bty='n',inset=c(0, -0.15), cex=cx) 

dev.off()

postscript("Fig3f_r0.7lambda1.7_Stoch.eps")
par(mai=c(2.5,0.9,0.4,0.8)+0.1, xpd=TRUE)

plot(gammaVector, gammaVector, type='n', lty=2, lwd=2, col=colA, 
     xlim=c(0,1), ylim=c(0,1), xlab=expression(paste("search-effort (", gamma, ")", sep='')), 
     ylab="proportion of area searched (n/N)",
     xaxs='i', yaxs='i',  xaxt="n", yaxt="n", cex.lab=1.4)
labelXY<-seq(0,1,0.2)
labelXY[1]<-0
labelXY[length(labelXY)]<-1
axis(side=1, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)
axis(side=2, at=seq(0,1,0.2), labels=labelXY, cex.axis=1.2)

for(i in 1:length(plotSurveys_II))
{  
  n_det[i,]<-numSUSearched_rApprox(p, pof, 0.7, plotSurveys_II[i], N, gammaVector) #r=0.7
  detCurve<-cbind(gammaVector, n_det[i,])
  
  if(length(which(detCurve[,2]>N)!=0))
  {
    detCurve<-detCurve[-which(detCurve[,2]>N),] #remove rows where n >N
  }
  
  detCurve[,2]<-detCurve[,2]/N #change to prp
  lines(detCurve[,1], detCurve[,2], lty=i+1, lwd=2, col=colA)
}

stochData<-solveForDelta(1.7, p, pof, kappaVector) #lambda=1.7
stochData<-as.data.frame(stochData)
colnames(stochData)<-c("kappa", "gammaN")

for(i in 1:length(plotSurveys_II))
{  
  n_stoch[i,]<-N*stochData[which(stochData$kappa==plotSurveys_II[i]),]$gammaN/gammaVector
  stochCurve<-cbind(gammaVector, n_stoch[i,])
  
  if(length(which(stochCurve[,2]>N)!=0))
  {
    stochCurve<-stochCurve[-which(stochCurve[,2]>N),] #remove rows where n >N
  }
  
  stochCurve[,2]<-stochCurve[,2]/N #change to prp
  lines(stochCurve[,1], stochCurve[,2], lty=i+1, lwd=2, col=colB)
  
}

legend("topright", legend=legText_II, lty=1, lwd=2, bty='n',inset=c(0, -0.15), cex=cx) 

dev.off()

