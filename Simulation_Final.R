args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  BATCH <- as.numeric(args[1])  #folder number
  set.seed(BATCH)
} else {
  stop()
}

PARA.ITER <- expand.grid( type=c(1,2), # 1=linear, 2=cubic
                          s=c(100,200,500,1000),
                          cut1=20,
                          cut2=40,
                          iter=c(1:10000))

source("Simulation_Functions.R")

PLOT <- FALSE
Var.Y <- 10
Boot.Num <- 1000 # number of bootstrap samples

GAMMA <- 1
KNOT <- c(10,10,10)

absbias.com <- function(x, tau.m)  apply(x, 2, function(x) abs(mean(x - tau.m)))
sd.com <- function(x)  apply(x, 2, function(x) sd(x))
RMSE.com <- function(x, tau.m)  apply(x, 2, function(x) sqrt(mean((x - tau.m)^2)))

outcome.type <- PARA.ITER[BATCH,1] # outcome model
s <- PARA.ITER[BATCH,2] # sample size
cut1 <- PARA.ITER[BATCH,3]  # regular testing time
cut2 <- PARA.ITER[BATCH,4] # max testing time
iter <- PARA.ITER[BATCH,5]

File.NDE <- sprintf("Summary_NDE_Type%0.1d_N_%0.4d_C1_%0.4d_C2_%0.4d_Iter_%0.5d.csv",
                    outcome.type,s,cut1,cut2,iter)
File.NIE <- sprintf("Summary_NIE_Type%0.1d_N_%0.4d_C1_%0.4d_C2_%0.4d_Iter_%0.5d.csv",
                    outcome.type,s,cut1,cut2,iter)
File.TE <- sprintf("Summary_TE_Type%0.1d_N_%0.4d_C1_%0.4d_C2_%0.4d_Iter_%0.5d.csv",
                   outcome.type,s,cut1,cut2,iter)
File.Cover <- sprintf("Summary_Cover_Type%0.1d_N_%0.4d_C1_%0.4d_C2_%0.4d_Iter_%0.5d.csv",
                      outcome.type,s,cut1,cut2,iter)
File.HP <- sprintf("Summary_HP_Type%0.1d_N_%0.4d_C1_%0.4d_C2_%0.4d_Iter_%0.5d.csv",
                   outcome.type,s,cut1,cut2,iter)

# DGP ---

k0 <- 5
if(outcome.type==1){
  
  Yk0.threshold <- Y.Z1.B.Ft.Lin(k0,0,1,cut1)
  dat <- data_generation(N=s, 
                         k0=k0, 
                         Yk0.threshold = Yk0.threshold, 
                         CUT=cut1, 
                         CUT2=cut2, 
                         Var.Y=Var.Y,
                         type="Lin")
} else {
  
  Yk0.threshold <- Y.Z1.B.Ft.Cub(k0,0,1,cut1)
  dat <- data_generation(N=s, 
                         k0=k0, 
                         Yk0.threshold = Yk0.threshold, 
                         CUT=cut1, 
                         CUT2=cut2, 
                         Var.Y=Var.Y,
                         type="Cub")
}

Data.All <- dat$Data.All
Data.Obs <- dat$Data.Obs

covariate_names <- c("X1", "X2")

# True Effects and Oracle Estimates

Oracle <- aggregate(cbind(Y.Z1.T1,Y.Z1.T1.Obs,
                          Y.Z1.TZ,Y.Z1.TZ.Obs,
                          Y.Z1.T0,Y.Z1.T0.Obs,
                          Y.Z0.T0,Y.Z0.T0.Obs)~K,Data.All,"mean")

True <- aggregate(cbind(Y.Z1.T1.True,
                        Y.Z1.TZ.True,
                        Y.Z1.T0.True,
                        Y.Z0.T0.True)~K,dat$Data.All,"mean")

NDE.Oracle <- Oracle$Y.Z1.T0.Obs - Oracle$Y.Z0.T0.Obs
NIE.Oracle <- Oracle$Y.Z1.TZ.Obs - Oracle$Y.Z1.T0.Obs
TE.Oracle  <- Oracle$Y.Z1.TZ.Obs - Oracle$Y.Z0.T0.Obs

NDE.True <- True$Y.Z1.T0.True - True$Y.Z0.T0.True
NIE.True <- True$Y.Z1.TZ.True - True$Y.Z1.T0.True
TE.True  <- True$Y.Z1.TZ.True - True$Y.Z0.T0.True

Y.Z0.T0.True <- True$Y.Z0.T0.True
Y.Z1.T0.True <- True$Y.Z1.T0.True
Y.Z1.TZ.True <- True$Y.Z1.TZ.True
Y.Z1.T1.True <- True$Y.Z1.T1.True

GAM.FULL.RESULT <-  Analysis(Data=Data.Obs,
                             covariate = covariate_names,
                             k0=k0,
                             Yk0="Yk0",
                             CUT=cut1,
                             CUT2=cut2,
                             x0=NULL,
                             type="GAM",
                             knot=KNOT,
                             gamma=GAMMA)

# Data=Data.Obs
# covariate = covariate_names
# k0=k0
# Yk0="Yk0"
# CUT=cut1
# CUT2=cut2
# x0=NULL
# type="GAM"
# knot=KNOT
# gamma=GAMMA

RESULT.Point.Estimate.Raw <- GAM.FULL.RESULT$RESULT
GAM.KNOT <- GAM.FULL.RESULT$knot


RESULT.Point.Estimate.K <- aggregate(.~K,RESULT.Point.Estimate.Raw,"mean")

NDE.GAM <- (RESULT.Point.Estimate.K$Y.Z1.T0 - RESULT.Point.Estimate.K$Y.Z0.T0)
NIE.GAM <- (RESULT.Point.Estimate.K$Y.Z1.TZ - RESULT.Point.Estimate.K$Y.Z1.T0)
TE.GAM  <- (RESULT.Point.Estimate.K$Y.Z1.TZ - RESULT.Point.Estimate.K$Y.Z0.T0)

# Linear approach

RESULT.Point.Estimate.Raw_Linear <- Analysis(Data=Data.Obs,
                                             covariate = covariate_names,
                                             k0=k0, 
                                             Yk0="Yk0",
                                             CUT=cut1,
                                             CUT2=cut2,
                                             x0=NULL,
                                             type="Linear")$RESULT
RESULT.Point.Estimate.K_Linear <- 
  aggregate(.~K,RESULT.Point.Estimate.Raw_Linear,"mean")

NDE.Linear <- (RESULT.Point.Estimate.K_Linear$Y.Z1.T0 - RESULT.Point.Estimate.K_Linear$Y.Z0.T0)
NIE.Linear <- (RESULT.Point.Estimate.K_Linear$Y.Z1.TZ - RESULT.Point.Estimate.K_Linear$Y.Z1.T0)
TE.Linear  <- (RESULT.Point.Estimate.K_Linear$Y.Z1.TZ - RESULT.Point.Estimate.K_Linear$Y.Z0.T0)

# Cubic approach

RESULT.Point.Estimate.Raw_Cubic <- Analysis(Data=Data.Obs,
                                            covariate = covariate_names,
                                            k0=k0, 
                                            Yk0="Yk0",
                                            CUT=cut1,
                                            CUT2=cut2,
                                            x0=NULL,
                                            type="Cubic")$RESULT
RESULT.Point.Estimate.K_Cubic <- 
  aggregate(.~K,RESULT.Point.Estimate.Raw_Cubic,"mean")

NDE.Cubic <- (RESULT.Point.Estimate.K_Cubic$Y.Z1.T0 - RESULT.Point.Estimate.K_Cubic$Y.Z0.T0)
NIE.Cubic <- (RESULT.Point.Estimate.K_Cubic$Y.Z1.TZ - RESULT.Point.Estimate.K_Cubic$Y.Z1.T0)
TE.Cubic  <- (RESULT.Point.Estimate.K_Cubic$Y.Z1.TZ - RESULT.Point.Estimate.K_Cubic$Y.Z0.T0)

# Bias

Oracle.Bias <- cbind(NDE.Oracle - NDE.True,
                     NIE.Oracle - NIE.True,
                     TE.Oracle - TE.True)

Linear.Bias <- cbind(NDE.Linear - NDE.True,
                     NIE.Linear - NIE.True,
                     TE.Linear - TE.True)

Cubic.Bias <- cbind(NDE.Cubic - NDE.True,
                    NIE.Cubic - NIE.True,
                    TE.Cubic - TE.True)

GAM.Bias <- cbind(NDE.GAM - NDE.True,
                  NIE.GAM - NIE.True,
                  TE.GAM - TE.True)



# Bootstrap

RESULT.Y.Z0.T0.Boot <- 
  RESULT.Y.Z1.T0.Boot <- 
  RESULT.Y.Z1.TZ.Boot <- 
  RESULT.Y.Z1.T1.Boot <- matrix(0,cut2,Boot.Num)

for(bb in 1:Boot.Num){
  
  set.seed(bb)
  ID.Boot <- sample(1:s,s,replace=TRUE)
  ID.Boot.unique <- unique(ID.Boot)
  ID.Boot.count <- table(factor(ID.Boot, levels = unique(ID.Boot)))
  Data.Obs.Boot <- do.call(rbind, replicate(ID.Boot.count[[1]], Data.Obs[Data.Obs$ID==ID.Boot.unique[1],], simplify = FALSE))
  for(i in 2:length(ID.Boot.unique)){
    Data.Obs.Boot2 <- do.call(rbind, replicate(ID.Boot.count[[i]], Data.Obs[Data.Obs$ID==ID.Boot.unique[i],], simplify = FALSE))
    Data.Obs.Boot <- rbind(Data.Obs.Boot, Data.Obs.Boot2)
  }
  
  RESULT.Point.Estimate.Raw.Boot <-
    tryCatch({Analysis(Data=Data.Obs.Boot,
                       covariate = covariate_names,
                       k0=k0,
                       Yk0="Yk0",
                       CUT=cut1,
                       CUT2=cut2,
                       x0=NULL,
                       type="GAM",
                       knot=KNOT,
                       gamma=GAMMA)},
             error = function(e) {})
  
  if(!is.null(RESULT.Point.Estimate.Raw.Boot)){
    RESULT.Y.Z0.T0.Boot[,bb] <- RESULT.Point.Estimate.Raw.Boot$RESULT$Y.Z0.T0
    RESULT.Y.Z1.T0.Boot[,bb] <- RESULT.Point.Estimate.Raw.Boot$RESULT$Y.Z1.T0
    RESULT.Y.Z1.TZ.Boot[,bb] <- RESULT.Point.Estimate.Raw.Boot$RESULT$Y.Z1.TZ
    RESULT.Y.Z1.T1.Boot[,bb] <- RESULT.Point.Estimate.Raw.Boot$RESULT$Y.Z1.T1
  } else {
    RESULT.Y.Z0.T0.Boot[,bb] <- NA
    RESULT.Y.Z1.T0.Boot[,bb] <- NA
    RESULT.Y.Z1.TZ.Boot[,bb] <- NA
    RESULT.Y.Z1.T1.Boot[,bb] <- NA
  }
  if(bb%%10==0){print(bb)}
  
}

NA.C <- which(apply(is.na(rbind(RESULT.Y.Z1.T0.Boot,RESULT.Y.Z0.T0.Boot,RESULT.Y.Z1.TZ.Boot)),2,sum)>0)
NNA.C <- setdiff(1:Boot.Num,NA.C)

NDE.Boot <- NIE.Boot <- TE.Boot <- 
  Y.Z0.T0.Boot <- Y.Z1.T0.Boot <- Y.Z1.T1.Boot <- Y.Z1.TZ.Boot <- 
  matrix(0,cut2,Boot.Num)

NDE.Boot[,NNA.C]     <- (RESULT.Y.Z1.T0.Boot - RESULT.Y.Z0.T0.Boot)[,NNA.C]
NIE.Boot[,NNA.C]     <- (RESULT.Y.Z1.TZ.Boot - RESULT.Y.Z1.T0.Boot)[,NNA.C]
TE.Boot[,NNA.C]      <- (RESULT.Y.Z1.TZ.Boot - RESULT.Y.Z0.T0.Boot)[,NNA.C]
Y.Z0.T0.Boot[,NNA.C] <- RESULT.Y.Z0.T0.Boot[,NNA.C]
Y.Z1.T0.Boot[,NNA.C] <- RESULT.Y.Z1.T0.Boot[,NNA.C]
Y.Z1.T1.Boot[,NNA.C] <- RESULT.Y.Z1.T1.Boot[,NNA.C]
Y.Z1.TZ.Boot[,NNA.C] <- RESULT.Y.Z1.TZ.Boot[,NNA.C]

NDE.Boot[,NA.C]     <- sample(c(-Inf,Inf),length(NA.C),replace=TRUE)
NIE.Boot[,NA.C]     <- sample(c(-Inf,Inf),length(NA.C),replace=TRUE) 
TE.Boot[,NA.C]      <- sample(c(-Inf,Inf),length(NA.C),replace=TRUE)
Y.Z0.T0.Boot[,NA.C] <- sample(c(-Inf,Inf),length(NA.C),replace=TRUE)
Y.Z1.T0.Boot[,NA.C] <- sample(c(-Inf,Inf),length(NA.C),replace=TRUE)
Y.Z1.T1.Boot[,NA.C] <- sample(c(-Inf,Inf),length(NA.C),replace=TRUE)
Y.Z1.TZ.Boot[,NA.C] <- sample(c(-Inf,Inf),length(NA.C),replace=TRUE)

NDE.CI <- apply(NDE.Boot,1,function(v){quantile(v,c(0.025,0.975))})
NIE.CI <- apply(NIE.Boot,1,function(v){quantile(v,c(0.025,0.975))})
TE.CI <- apply(TE.Boot,1,function(v){quantile(v,c(0.025,0.975))})
Y.Z0.T0.CI <- apply(Y.Z0.T0.Boot,1,function(v){quantile(v,c(0.025,0.975))})
Y.Z1.T0.CI <- apply(Y.Z1.T0.Boot,1,function(v){quantile(v,c(0.025,0.975))})
Y.Z1.T1.CI <- apply(Y.Z1.T1.Boot,1,function(v){quantile(v,c(0.025,0.975))})
Y.Z1.TZ.CI <- apply(Y.Z1.TZ.Boot,1,function(v){quantile(v,c(0.025,0.975))})

Cover <- cbind( as.numeric(NDE.CI[1,] <= NDE.True & NDE.True <= NDE.CI[2,]),
                as.numeric(NIE.CI[1,] <= NIE.True & NIE.True <= NIE.CI[2,]),
                as.numeric(TE.CI[1,] <= TE.True & TE.True <= TE.CI[2,]),
                as.numeric(Y.Z0.T0.CI[1,] <= Y.Z0.T0.True & Y.Z0.T0.True <= Y.Z0.T0.CI[2,]),
                as.numeric(Y.Z1.T0.CI[1,] <= Y.Z1.T0.True & Y.Z1.T0.True <= Y.Z1.T0.CI[2,]),
                as.numeric(Y.Z1.T1.CI[1,] <= Y.Z1.T1.True & Y.Z1.T1.True <= Y.Z1.T1.CI[2,]),
                as.numeric(Y.Z1.TZ.CI[1,] <= Y.Z1.TZ.True & Y.Z1.TZ.True <= Y.Z1.TZ.CI[2,]))

# Summary

NDE.Summary <- c( Oracle.Bias[,1],
                  Linear.Bias[,1],
                  Cubic.Bias[,1],
                  GAM.Bias[,1],
                  Cover[,1])
NDE.Summary <- matrix(NDE.Summary,1,length(NDE.Summary))

NIE.Summary <- c( Oracle.Bias[,2],
                  Linear.Bias[,2],
                  Cubic.Bias[,2],
                  GAM.Bias[,2],
                  Cover[,2])
NIE.Summary <- matrix(NIE.Summary,1,length(NIE.Summary))

TE.Summary <- c( Oracle.Bias[,3],
                 Linear.Bias[,3],
                 Cubic.Bias[,3],
                 GAM.Bias[,3],
                 Cover[,3])
TE.Summary <- matrix(TE.Summary,1,length(TE.Summary))

colnames(NDE.Summary) <- 
  colnames(NIE.Summary) <- 
  colnames(TE.Summary) <- c( sprintf("Bias_Oracle_%0.4d",1:cut2),
                             sprintf("Bias_Linear_%0.4d",1:cut2),
                             sprintf("Bias_Cubic_%0.4d",1:cut2),
                             sprintf("Bias_GAM_%0.4d",1:cut2),
                             sprintf("Cover_GAM_%0.4d",1:cut2))

Cover.Basis <- c(Cover[,1],Cover[,2],Cover[,3],
                 Cover[,4],Cover[,5],Cover[,6],Cover[,7])
Cover.Basis <- matrix(Cover.Basis,1,length(Cover.Basis))
colnames(Cover.Basis) <- c( sprintf("Cover_NDE_%0.4d",1:cut2),
                            sprintf("Cover_NIE_%0.4d",1:cut2),
                            sprintf("Cover_TE_%0.4d",1:cut2),
                            sprintf("Cover_Y.Z0.T0_%0.4d",1:cut2),
                            sprintf("Cover_Y.Z1.T0_%0.4d",1:cut2),
                            sprintf("Cover_Y.Z1.T1_%0.4d",1:cut2),
                            sprintf("Cover_Y.Z1.TZ_%0.4d",1:cut2) )



VT <- (k0+1):cut2

if(PLOT){
  
  layout(matrix(1:16,2,8,byrow=T))
  par(mar=c(3,3,1,0.5))
  plot(VT,Y.Z0.T0.True[VT],type='l',col=1,lwd=2,ylim=c(0,max(True[,2:5])))
  points(VT,Y.Z1.T0.True[VT],type='l',col=2,lwd=2)
  points(VT,Y.Z1.T1.True[VT],type='l',col=3,lwd=2)
  points(VT,Y.Z1.TZ.True[VT],type='l',col=4,lwd=2)
  
  points(VT,Oracle$Y.Z0.T0.Obs[VT],type='l',col=1,lwd=1)
  points(VT,Oracle$Y.Z1.T0.Obs[VT],type='l',col=2,lwd=1)
  points(VT,Oracle$Y.Z1.T1.Obs[VT],type='l',col=3,lwd=1)
  points(VT,Oracle$Y.Z1.TZ.Obs[VT],type='l',col=4,lwd=1)
  
  plot(VT,Y.Z0.T0.True[VT],type='l',col=1, lwd=3,ylim=c(0,max(True[,5])+5))
  points(VT,RESULT.Point.Estimate.K$Y.Z0.T0[VT],type='l',col=1,lwd=1)
  
  plot(VT,Y.Z1.T0.True[VT],type='l',col=2, lwd=3,ylim=c(0,max(True[,4])+5))
  points(VT,RESULT.Point.Estimate.K$Y.Z1.T0[VT],type='l',col=2,lwd=1)
  
  plot(VT,Y.Z1.TZ.True[VT],type='l',col=4, lwd=3,ylim=c(0,max(True[,3])+5))
  points(VT,RESULT.Point.Estimate.K$Y.Z1.TZ[VT],type='l',col=4,lwd=1)
  
  plot(VT,Y.Z1.T1.True[VT],type='l',col=3, lwd=2,ylim=c(0,max(True[,2])+5))
  points(VT,RESULT.Point.Estimate.K$Y.Z1.T1[VT],type='l',col=3,lwd=1)
  
  plot(VT,NDE.True[VT],type='l',col=4, lwd=2)
  points(VT,NDE.GAM[VT],type='l',col=2,lwd=2)
  points(VT,NDE.CI[1,VT],type='l',col=1,lwd=1,lty=2)
  points(VT,NDE.CI[2,VT],type='l',col=1,lwd=1,lty=2)
  
  plot(VT,NIE.True[VT],type='l',col=4, lwd=2)
  points(VT,NIE.GAM[VT],type='l',col=2,lwd=2)
  points(VT,NIE.CI[1,VT],type='l',col=1,lwd=1,lty=2)
  points(VT,NIE.CI[2,VT],type='l',col=1,lwd=1,lty=2)
  
  
  plot(VT,TE.True[VT],type='l',col=4, lwd=2)
  points(VT,TE.GAM[VT],type='l',col=2,lwd=2)
  points(VT,TE.CI[1,VT],type='l',col=1,lwd=1,lty=2)
  points(VT,TE.CI[2,VT],type='l',col=1,lwd=1,lty=2)
  
  ###
  
  plot.new()
  
  plot(VT,(Y.Z0.T0.True -
             RESULT.Point.Estimate.K$Y.Z0.T0)[VT],type='l',col=1, lwd=2)
  
  plot(VT,(Y.Z1.T0.True -
             RESULT.Point.Estimate.K$Y.Z1.T0)[VT],type='l',col=2, lwd=2)
  
  plot(VT,(Y.Z1.TZ.True -
             RESULT.Point.Estimate.K$Y.Z1.TZ)[VT],type='l',col=4, lwd=2)
  
  plot(VT,(Y.Z1.T1.True -
             RESULT.Point.Estimate.K$Y.Z1.T1)[VT],type='l',col=3,lwd=2)
  
  plot(VT,(NDE.True - NDE.GAM)[VT],type='l',col=4, lwd=2, 
       ylim=range(c((NDE.True - NDE.CI[2,])[VT],
                    (NDE.True - NDE.CI[1,])[VT[length(VT):1]])))
  
  text(cut2,(NDE.True - NDE.CI[2,])[cut2],
       round(mean((NDE.True >= NDE.CI[1,] & NDE.True  <= NDE.CI[2,])[VT]),3),pos=2)
  
  polygon(cbind(c(VT,VT[length(VT):1]),
                c((NDE.True - NDE.CI[1,])[VT],
                  (NDE.True - NDE.CI[2,])[VT[length(VT):1]])),
          col=rgb(0,0,0,0.25),
          border=F)
  abline(h=0,col=1)
  
  plot(VT,(NIE.True - NIE.GAM)[VT],type='l',col=4, lwd=2, 
       ylim=range(c(NIE.True - NIE.CI[2,],
                    NIE.True - NIE.CI[1,])))
  
  text(cut2,(NIE.True - NIE.CI[2,])[cut2],
       round(mean((NIE.True >= NIE.CI[1,] & NIE.True  <= NIE.CI[2,])[VT]),3),pos=2)
  
  polygon(cbind(c(VT,VT[length(VT):1]),
                c((NIE.True - NIE.CI[1,])[VT],
                  (NIE.True - NIE.CI[2,])[VT[length(VT):1]])),
          col=rgb(0,0,0,0.25),
          border=F)
  abline(h=0,col=1)
  
  
  plot(VT,(TE.True - TE.GAM)[VT],type='l',col=4, lwd=2, 
       ylim=range(c(TE.True - TE.CI[2,],
                    TE.True - TE.CI[1,])))
  
  text(cut2,(TE.True - TE.CI[2,])[cut2],
       round(mean((TE.True >= TE.CI[1,] & TE.True  <= TE.CI[2,])[VT]),3),pos=2)
  
  polygon(cbind(c(VT,VT[length(VT):1]),
                c((TE.True - TE.CI[1,])[VT],
                  (TE.True - TE.CI[2,])[VT[length(VT):1]])),
          col=rgb(0,0,0,0.25),
          border=F)
  abline(h=0,col=1)
  
}

write.csv(NDE.Summary,File.NDE,row.names=F)
write.csv(NIE.Summary,File.NIE,row.names=F)
write.csv(TE.Summary,File.TE,row.names=F)
write.csv(Cover.Basis,File.Cover,row.names=F)
 