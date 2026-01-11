library(mgcv)
library(nnet)
library(tidyverse)

# ::: Analysis: GAM Model ::: ####

Analysis <- function(Data, 
                     k0=5,
                     Yk0="Yk0",
                     CUT=20, 
                     CUT2=40, 
                     covariate, 
                     x0=NULL,
                     type="GAM",
                     knot=NULL, 
                     gamma=1){
  
  # Data <- data_extend(Data)
  
  CN <- colnames(Data)
  CN.X <- CN[substr(CN,1,1)=="X"]
  CN.Yk0X <- c(Yk0,CN[substr(CN,1,1)=="X"])
  
  Data <- Data[order(Data[,CN.X[1]],Data[,CN.X[2]]),]
  Data <- Data[Data$K>k0,]                            # Throw away the beginning of the time series
  
  Data.Long <- Data
  Data.Short <- unique(Data[,c("ID","Z","T",CN.Yk0X)])
  Data.Yk0X.Unique <- as.data.frame(unique(Data.Short[,CN.Yk0X]))
  
  colnames(Data.Yk0X.Unique) <- CN.Yk0X
  Data.Yk0X.Unique$weight <- sapply(1:nrow(Data.Yk0X.Unique),function(ii){
    match <- Reduce('&',lapply(covariate, function(col){
      Data.Short[[col]]== Data.Yk0X.Unique[[col]][ii]
    }))
    mean(match)
  })
  
  Data.X.Unique <- as.data.frame(unique(Data.Short[,CN.X]))
  colnames(Data.X.Unique) <- CN.X
  Data.X.Unique$weight <- sapply(1:nrow(Data.X.Unique),function(ii){
    match <- Reduce('&',lapply(covariate, function(col){
      Data.Short[[col]]== Data.X.Unique[[col]][ii]
    }))
    mean(match)
  })
  
  if(type=="GAM"){
    
    Z0T0 <- which(Data.Long$Z==0 & Data.Long$T==0)
    Data.Long.Z0T0 <- Data.Long[Z0T0,]
    Z1T0 <- which(Data.Long$Z==1 & Data.Long$T==0)
    Data.Long.Z1T0 <- Data.Long[Z1T0,]
    Z1T1 <- which(Data.Long$Z==1 & Data.Long$T==1)
    Data.Long.Z1T1 <- Data.Long[Z1T1,]
    
    SELECT <- FALSE
    
    model.cut1 <- knot[1]
    model.cut2 <- knot[2]
    model.cut3 <- knot[3]
    
    covariates <- paste(covariate, collapse = "+")
    covariates_Yk0 <- paste(c(covariate,Yk0), collapse = "+")
    covariates_complex <- paste(covariate, collapse = "*")
    
    gam_intercept_model41 <- paste(sapply(c(covariate,Yk0), function(x) paste(sprintf("s(K,by=%s,k=%s)",x,model.cut1),sep = "")),collapse = "+")
    gam_intercept_model42 <- paste(sapply(c(covariate,Yk0), function(x) paste(sprintf("s(K,by=%s,k=%s)",x,model.cut2),sep = "")),collapse = "+")
    gam_intercept_model43 <- paste(sapply(c(covariate,Yk0), function(x) paste(sprintf("s(K,by=%s,k=%s)",x,model.cut3),sep = "")),collapse = "+")

    gam_model_1 <- as.formula(paste(sprintf("Y~s(K,k=%s)+",model.cut1),
                                    covariates_Yk0,"+",
                                    gam_intercept_model41))
    gam_model_2 <- as.formula(paste(sprintf("Y~s(K,k=%s)+",model.cut2),
                                    covariates_Yk0,"+",
                                    gam_intercept_model42))
    gam_model_3 <- as.formula(paste(sprintf("Y~s(K,k=%s)+",model.cut3),
                                    covariates_Yk0,"+",
                                    gam_intercept_model43))
    
    gam_model_short_1 <- as.formula(paste(sprintf("Y~1+s(K,k=%s)+",model.cut1), covariates_Yk0))
    gam_model_short_2 <- as.formula(paste(sprintf("Y~1+s(K,k=%s)+",model.cut2), covariates_Yk0))
    gam_model_short_3 <- as.formula(paste(sprintf("Y~1+s(K,k=%s)+",model.cut3), covariates_Yk0))
    gam_model_T <- as.formula(paste("T~",covariates_complex))
    gam_model_Yk0 <- as.formula(paste("Yk0~",covariates_complex))
    
    GAM.Y.Z0.T0 <- try(gam(gam_model_1, data=Data.Long.Z0T0,select=SELECT,gamma=gamma,method="REML"),silent=TRUE)
    GAM.Y.Z1.T0 <- try(gam(gam_model_2, data=Data.Long.Z1T0,select=SELECT,gamma=gamma,method="REML"),silent=TRUE)
    GAM.Y.Z1.T1 <- try(gam(gam_model_3, data=Data.Long.Z1T1,select=SELECT,gamma=gamma,method="REML"),silent=TRUE)

    if(class(GAM.Y.Z0.T0)[1]=="try-error"){
    GAM.Y.Z0.T0 <-
      try({gam(gam_model_short_1, data=Data.Long.Z0T0,select=SELECT,gamma=gamma,method="REML")},silent=TRUE)
    }
    if(class(GAM.Y.Z1.T0)[1]=="try-error"){
    GAM.Y.Z1.T0 <-
      try({gam(gam_model_short_2, data=Data.Long.Z1T0,select=SELECT,gamma=gamma,method="REML")},silent=TRUE)
    }
    if(class(GAM.Y.Z1.T1)[1]=="try-error"){
    GAM.Y.Z1.T1 <-
      try({gam(gam_model_short_3, data=Data.Long.Z1T1,select=SELECT,gamma=gamma,method="REML")},silent=TRUE)
    }
    
    ProbYk0_T1Z1 <- gam(gam_model_Yk0, 
                        data=Data.Short[Data.Short$Z==1 & Data.Short$T==1,], 
                        family="binomial",select=FALSE,gamma=gamma,method="REML")
    ProbYk0_T0Z1 <- gam(gam_model_Yk0, 
                        data=Data.Short[Data.Short$Z==1 & Data.Short$T==0,], 
                        family="binomial",select=FALSE,gamma=gamma,method="REML")
    ProbYk0_Z0 <- gam(gam_model_Yk0, data=Data.Short[Data.Short$Z==0,], 
                      family="binomial",select=FALSE,gamma=gamma,method="REML")
    ProbT <- gam(gam_model_T, data=Data.Short[Data.Short$Z==1,], 
                 family="binomial",select=FALSE,gamma=gamma,method="REML")
    
  } else if(type=="Linear") {
    
    model.cut1 <- model.cut2 <- model.cut3 <- NULL
    covariates <- paste(covariate, collapse = "+")
    covariates_Yk0 <- paste(c(covariate,Yk0), collapse = "+")
    covariates_complex <- paste(covariate, collapse = "*") 
    
    lm_interaction_model <- as.formula(paste("Y~K+", 
                                             covariates_Yk0,"+",
                                             paste(sapply(c(covariate,Yk0), function(x) paste("K:",x,"",sep = "")),
                                                   collapse = "+")))
    glm_Yk0_model <- as.formula(paste("Yk0~",covariates_complex)) 
    glm_model <- as.formula(paste("T~",covariates_complex))
    
    GAM.Y.Z0.T0 <- lm(lm_interaction_model, data=Data.Long[Data.Long$Z==0 & Data.Long$T==0,])
    GAM.Y.Z1.T0 <- lm(lm_interaction_model, data=Data.Long[Data.Long$Z==1 & Data.Long$T==0,])
    GAM.Y.Z1.T1 <- lm(lm_interaction_model, data=Data.Long[Data.Long$Z==1 & Data.Long$T==1,])
    
    # Y ~ K + X1 + X2 + Yk0 + K:X1 + K:X2 + K:Yk0
    
    ProbYk0_T1Z1 <- glm(glm_Yk0_model, 
                        data=Data.Short[Data.Short$Z==1 & Data.Short$T==1,], 
                        family="binomial")
    ProbYk0_T0Z1 <- glm(glm_Yk0_model, 
                        data=Data.Short[Data.Short$Z==1 & Data.Short$T==0,], 
                        family="binomial")
    ProbYk0_Z0 <- glm(glm_Yk0_model, 
                      data=Data.Short[Data.Short$Z==0,], 
                      family="binomial")
    ProbT <- glm(glm_model, data=Data.Short[Data.Short$Z==1,], family="binomial")
    
  } else if(type=="Cubic") {
    
    model.cut1 <- model.cut2 <- model.cut3 <- NULL
    covariates <- paste(covariate, collapse = "+")
    covariates_Yk0 <- paste(c(covariate,Yk0), collapse = "+")
    covariates_complex <- paste(covariate, collapse = "*") 
    
    lm_interaction_model <- as.formula(paste("Y~K+I(K^2)+I(K^3)+", covariates_Yk0,"+",
                                             paste(sapply(c(covariate,Yk0), function(x) paste("K:",x,"",sep = "")),
                                                   collapse = "+")))
    glm_Yk0_model <- as.formula(paste("Yk0~",covariates_complex)) 
    glm_model <- as.formula(paste("T~",covariates_complex))
    
    GAM.Y.Z0.T0 <- lm(lm_interaction_model, data=Data.Long[Data.Long$Z==0 & Data.Long$T==0,])
    GAM.Y.Z1.T0 <- lm(lm_interaction_model, data=Data.Long[Data.Long$Z==1 & Data.Long$T==0,])
    GAM.Y.Z1.T1 <- lm(lm_interaction_model, data=Data.Long[Data.Long$Z==1 & Data.Long$T==1,])
    
    # K + I(K^2) + I(K^3) + X1 + X2 + Yk0 + K:X1 + K:X2 + K:Yk0
    
    ProbYk0_T1Z1 <- glm(glm_Yk0_model, 
                        data=Data.Short[Data.Short$Z==1 & Data.Short$T==1,], 
                        family="binomial")
    ProbYk0_T0Z1 <- glm(glm_Yk0_model, 
                        data=Data.Short[Data.Short$Z==1 & Data.Short$T==0,], 
                        family="binomial")
    ProbYk0_Z0 <- glm(glm_Yk0_model, 
                      data=Data.Short[Data.Short$Z==0,], 
                      family="binomial")
    ProbT <- glm(glm_model, data=Data.Short[Data.Short$Z==1,], family="binomial")
    
    
  }
  
  if(is.null(x0)){
    
    Newdata.KX <- data.frame(matrix(0,CUT2*nrow(Data.X.Unique),length(covariate)+2)) ####
    colnames(Newdata.KX) <- c("K",covariate,"weight") ####
    
    Newdata.KX$K <- 1:CUT2
    for(nr in 1:nrow(Data.X.Unique)){ #######
      Newdata.KX[((nr-1)*CUT2+1):(nr*CUT2),2:(length(covariate)+1)] <- 
        matrix(rep(as.numeric(Data.X.Unique[nr,1:length(covariate)]), CUT2),
               nrow=CUT2, ncol=length(covariate), byrow=TRUE)
    }
    
    if(length(covariate)==1){
      Newdata.KX$weight <- 
        sapply(1:nrow(Newdata.KX),
               function(ii){
                 Data.X.Unique$weight[which(Newdata.KX$X1[ii]==Data.X.Unique$X1)]
               })*nrow(Data.X.Unique)
    } else{
      Newdata.KX[seq(1,(nrow(Newdata.KX)-CUT2+1),CUT2),]$weight <- 
        sapply(seq(1,(nrow(Newdata.KX)-CUT2+1),CUT2), function(ii){
          Data.X.Unique$weight[apply(Data.X.Unique,1,function(row) all(row==Newdata.KX[ii, CN.X]))]
        }) * nrow(Data.X.Unique)
      Newdata.KX <- Newdata.KX %>% 
        mutate(weight = ifelse(weight == 0, NA, weight)) %>% 
        fill(weight, .direction = "down")
      
      
    }
    
  } else {
    list <- list(K = 1:CUT2)
    for(i in seq_along(x0)){
      list[[paste0("X",i)]] <- x0[i]
    }
    Newdata.KX <- do.call(expand.grid, list)
    Newdata.KX$weight <- 1
  }
  
  Newdata.KX.CUT <- Newdata.KX
  Newdata.KX.CUT[Newdata.KX.CUT$K > CUT, ]$K <- CUT
  
  Newdata.KXYk0 <- rbind( cbind(Newdata.KX,0),
                          cbind(Newdata.KX,0) )
  Newdata.KXYk0[nrow(Newdata.KX)+1:nrow(Newdata.KX),ncol(Newdata.KXYk0)] <- 1
  colnames(Newdata.KXYk0) <- c(colnames(Newdata.KX),Yk0)
  
  Newdata.KXYk0.CUT <- Newdata.KXYk0
  Newdata.KXYk0.CUT[Newdata.KXYk0.CUT$K > CUT,]$K <- CUT
  
  Newdata.KXYk0.T1Z1 <- Newdata.KXYk0.T0Z1 <- Newdata.KXYk0.Z0 <- Newdata.KXYk0
  Newdata.KXYk0.T1Z1.CUT <- Newdata.KXYk0.T0Z1.CUT <- Newdata.KXYk0.Z0.CUT <- Newdata.KXYk0.CUT
  
  ## Pr[Yk0=1|Z=z,T=t,X]
  ProbYk01_1_T1Z1 <- predict(ProbYk0_T1Z1,newdata=Newdata.KXYk0.T1Z1,type="response")
  ProbYk01_0_T1Z1 <- 1-ProbYk01_1_T1Z1
  ProbYk01_1_T0Z1 <- predict(ProbYk0_T0Z1,newdata=Newdata.KXYk0.T0Z1,type="response")
  ProbYk01_0_T0Z1 <- 1-ProbYk01_1_T0Z1
  ProbYk01_1_Z0 <- predict(ProbYk0_Z0,newdata=Newdata.KXYk0.Z0,type="response")
  ProbYk01_0_Z0 <- 1-ProbYk01_1_Z0
  
  ProbYk01_y_T1Z1 <- (ProbYk01_1_T1Z1*Newdata.KXYk0.T1Z1$Yk0 + ProbYk01_0_T1Z1*(1-Newdata.KXYk0.T1Z1$Yk0))*2
  ProbYk01_y_T0Z1 <- (ProbYk01_1_T0Z1*Newdata.KXYk0.T0Z1$Yk0 + ProbYk01_0_T0Z1*(1-Newdata.KXYk0.T0Z1$Yk0))*2
  ProbYk01_y_Z0 <- (ProbYk01_1_Z0*Newdata.KXYk0.Z0$Yk0 + ProbYk01_0_Z0*(1-Newdata.KXYk0.Z0$Yk0))*2
  
  Newdata.KXYk0.T1Z1$weight.y <- Newdata.KXYk0.T1Z1$weight*ProbYk01_y_T1Z1
  Newdata.KXYk0.T0Z1$weight.y <- Newdata.KXYk0.T0Z1$weight*ProbYk01_y_T0Z1
  Newdata.KXYk0.Z0$weight.y <- Newdata.KXYk0.Z0$weight*ProbYk01_y_Z0
  
  ## Pr[T=1|Z=1,X]
  ProbT1 <- predict(ProbT,newdata=Newdata.KXYk0.T1Z1,type="response")
  
  # ProbT1 <- rep(c(0.7,0.6,0.6,0.5),each=40)
  # Newdata.KX$weight <- 1
  
  
  ## E[Y|Z=0,T=0,H,X]
  
  muY.Z0.T0.GAM.Mat <- data.frame(K=Newdata.KXYk0.Z0$K,
                                  Yk0=Newdata.KXYk0.Z0$Yk0,
                                  Pred=predict(GAM.Y.Z0.T0, 
                                               newdata=Newdata.KXYk0.Z0.CUT, 
                                               type="response"))
  muY.Z0.T0.GAM.Mat$Pred.Weight <- muY.Z0.T0.GAM.Mat$Pred*Newdata.KXYk0.Z0$weight.y 
  
  muY.Z0.T0.GAM <- aggregate(Pred.Weight~K,muY.Z0.T0.GAM.Mat,"mean")$Pred.Weight
  
  ## E[Y|Z=1,T=0,H,X]
  
  muY.Z1.T0.GAM.Mat <- data.frame(K=Newdata.KXYk0.T0Z1$K,
                                  Pred=predict(GAM.Y.Z1.T0,
                                               newdata=Newdata.KXYk0.T0Z1.CUT, 
                                               type="response"))
  
  muY.Z1.T0.GAM.Mat$Pred.Weight <- muY.Z1.T0.GAM.Mat$Pred*Newdata.KXYk0.T0Z1$weight.y
  
  muY.Z1.T0.GAM <- aggregate(Pred.Weight~K,muY.Z1.T0.GAM.Mat,"mean")$Pred.Weight
  
  ## E[Y|Z=1,T=1,H,X]
  
  muY.Z1.T1.GAM.Mat <- data.frame(K=Newdata.KXYk0.T1Z1$K,
                                  Pred=predict(GAM.Y.Z1.T1,
                                               newdata=Newdata.KXYk0.T1Z1,
                                               type="response"))
  
  muY.Z1.T1.GAM.Mat$Pred.Weight <- muY.Z1.T1.GAM.Mat$Pred*Newdata.KXYk0.T1Z1$weight.y
  
  muY.Z1.T1.GAM <- aggregate(Pred.Weight~K,muY.Z1.T1.GAM.Mat,"mean")$Pred.Weight
  
  
  muY.Z1.TZ.GAM.Mat <- muY.Z1.T1.GAM.Mat
  muY.Z1.TZ.GAM.Mat$Pred.Weight <- 
    (muY.Z1.T1.GAM.Mat$Pred*Newdata.KXYk0.T1Z1$weight.y*ProbT1 + 
       muY.Z1.T0.GAM.Mat$Pred*Newdata.KXYk0.T0Z1$weight.y*(1-ProbT1))
  # muY.Z1.TZ.GAM.Mat$Pred.Weight <- 
  #   muY.Z1.TZ.GAM.Mat$Pred.TYk0*Newdata.KXYk0.Z1$weight
  muY.Z1.TZ.GAM <- aggregate(Pred.Weight~K,muY.Z1.TZ.GAM.Mat,"mean")$Pred.Weight
  
  output <- list()
  RESULT <- data.frame(cbind(muY.Z0.T0.GAM,
                             muY.Z1.T0.GAM,
                             muY.Z1.TZ.GAM,
                             muY.Z1.T1.GAM,
                             ProbT1[1:CUT2],
                             Newdata.KX$K[1:CUT2],
                             Newdata.KX[1:CUT2,CN.X]))
  colnames(RESULT)<- c(
    "Y.Z0.T0",
    "Y.Z1.T0",
    "Y.Z1.TZ",
    "Y.Z1.T1",
    "ProbT1",
    "K",
    CN.X)
  output$RESULT <- RESULT
  output$knot <- c(model.cut1,
                   model.cut2,
                   model.cut3) 
  return(output)
}


# ::: Data Generating Model ::: ####

data_generation <- function(N=1000, 
                            k0=5,
                            Yk0.threshold=10,
                            CUT=20, 
                            CUT2=40, 
                            Var.Y=10, 
                            noC_Zassign=FALSE,
                            type="Lin") {
  
  if(type=="Lin"){
    #############################
    # Potential Score
    #############################
    
    Y.Z0.B.Ft <- function(k,h,x){ Y.Z0.B.Ft.Lin(k,h,x,CUT) }
    Y.Z1.B.Ft <- function(k,h,x){ Y.Z1.B.Ft.Lin(k,h,x,CUT) }
    
  } else {
    #############################
    # Potential Score
    #############################
    
    Y.Z0.B.Ft <- function(k,h,x){ Y.Z0.B.Ft.Cub(k,h,x,CUT) }
    Y.Z1.B.Ft <- function(k,h,x){ Y.Z1.B.Ft.Cub(k,h,x,CUT) }
    
  }
  
  
  
  
  
  #############################
  # Potential Horizon
  #############################
  
  H0.Ft <- function(y,x){                                                              # expected potential outcome for whom used A=0 and H<=30 (i.e., non-users)
    if(y==0){  # outcome=low -> use more time
      if(x==0){
        sample(c(10,15,20),1,prob=c(1,1,3))
      } else if (x==1) {
        sample(c(10,15,20),1,prob=c(3,3,4))
      } else if (x==2) {
        sample(c(10,15,20),1,prob=c(1,1,1))
      }
    } else {  # outcome=high -> use less time
      if(x==0){
        sample(c(10,15,20),1,prob=c(1,1,2))
      } else if (x==1) {
        sample(c(10,15,20),1,prob=c(1,1,1))
      } else if (x==2) {
        sample(c(10,15,20),1,prob=c(1,1,1))
      }
    }
    
  }
  
  p00 <- c(1,1,3); p01 <- c(3,3,4); p02 <- c(1,1,1)
  p10 <- c(1,1,2); p11 <- c(1,1,1); p12 <- c(1,1,1)
  
  p00 <- p00/sum(p00)
  p01 <- p01/sum(p01)
  p02 <- p02/sum(p02)
  p10 <- p10/sum(p10)
  p11 <- p11/sum(p11)
  p12 <- p12/sum(p12)
  
  H1.Ft <- function(y,x){                                                              # expected potential outcome for whom used T>30 (i.e., users) 
    if(y==0){  # outcome=low -> use more time
      if(x==0){
        sample(c(20,30,40),1,prob=p00)
      } else if (x==1) {
        sample(c(20,30,40),1,prob=p01)
      } else if (x==2) {
        sample(c(20,30,40),1,prob=p02)
      }
    } else {
      if(x==0){
        sample(c(20,30,40),1,prob=p10)
      } else if (x==1) {
        sample(c(20,30,40),1,prob=p11)
      } else if (x==2) {
        sample(c(20,30,40),1,prob=p12)
      } 
    }
    
  }
  
  
  
  ################################################################################
  # Data for Estimation
  ################################################################################
  
  X1 <- rbinom(N,1,0.5)                                                                # reference value for X
  X2 <- rbinom(N,1,0.5) 
  
  Z <- rbinom(N,1,exp((X1+X2-1)*0.2)/(1+exp((X1+X2-1)*0.2))) 
  
  #############################
  # Outcome
  #############################
  
  XX <- c(0,1,1,2)
  
  prob.of.Yk0.1 <-
    (1-pnorm(Yk0.threshold,mean=Y.Z1.B.Ft(k0,0,0),sd=Var.Y))*0.25 +
    (1-pnorm(Yk0.threshold,mean=Y.Z1.B.Ft(k0,0,1),sd=Var.Y))*0.5 +
    (1-pnorm(Yk0.threshold,mean=Y.Z1.B.Ft(k0,0,2),sd=Var.Y))*0.25
  
  prob.of.Yk0.0 <- 1-prob.of.Yk0.1
  
  P0.Vec <- c( (p00[1])*prob.of.Yk0.0 + (p10[1])*prob.of.Yk0.1,
               (p01[1])*prob.of.Yk0.0 + (p11[1])*prob.of.Yk0.1,
               (p02[1])*prob.of.Yk0.0 + (p12[1])*prob.of.Yk0.1)
  
  P0.Vec.F <- P0.Vec[c(1,2,2,3)]
  
  # P0.Vec <- c(0.2,0.3,0.3,1/3)
  
  Y.Z1.T1.Mean.Long <- Y.Z1.TZ.Mean.Long <- Y.Z1.T0.Mean.Long <- Y.Z0.T0.Mean.Long <- rep(0,CUT2*N)
  for(ii in 1:N){
    for(jj in 1:CUT2){
      Y.Z0.T0.Mean.Long[(ii-1)*CUT2+jj] <- Y.Z0.B.Ft(jj,0,X1[ii]+X2[ii])
      Y.Z1.T0.Mean.Long[(ii-1)*CUT2+jj] <- Y.Z1.B.Ft(jj,0,X1[ii]+X2[ii])
      Y.Z1.T1.Mean.Long[(ii-1)*CUT2+jj] <- Y.Z1.B.Ft(jj,CUT2,X1[ii]+X2[ii])
      Y.Z1.TZ.Mean.Long[(ii-1)*CUT2+jj] <- 
        Y.Z1.T1.Mean.Long[(ii-1)*CUT2+jj]*(1-P0.Vec[X1[ii]+X2[ii]+1])+
        Y.Z1.T0.Mean.Long[(ii-1)*CUT2+jj]*(P0.Vec[X1[ii]+X2[ii]+1])
    }
  }
  
  
  Y.Z1.T1.True.Mat <- Y.Z1.TZ.True.Mat <- Y.Z1.T0.True.Mat <- Y.Z0.T0.True.Mat <- matrix(0,CUT2,4)
  for(jj in 1:CUT2){ 
    for(tt in 1:4){
      Y.Z0.T0.True.Mat[jj,tt] <- Y.Z0.B.Ft(jj,0,XX[tt])
      Y.Z1.T0.True.Mat[jj,tt] <- Y.Z1.B.Ft(jj,0,XX[tt])
      Y.Z1.T1.True.Mat[jj,tt] <- Y.Z1.B.Ft(jj,CUT2,XX[tt])
      Y.Z1.TZ.True.Mat[jj,tt] <- 
        (Y.Z1.T1.True.Mat[jj,tt])*(1-P0.Vec.F[tt]) + 
        (Y.Z1.T0.True.Mat[jj,tt])*(P0.Vec.F[tt])
    }
  }
  
  Y.Z1.T1.True <- apply(Y.Z1.T1.True.Mat,1,mean)
  Y.Z1.TZ.True <- apply(Y.Z1.TZ.True.Mat,1,mean)
  Y.Z1.T0.True <- apply(Y.Z1.T0.True.Mat,1,mean)
  Y.Z0.T0.True <- apply(Y.Z0.T0.True.Mat,1,mean)
  
  ID.Long <- rep(1:N,each=CUT2)
  Z.Long <- rep(Z,each=CUT2)
  Y.Z1.T1.Long <- Y.Z1.T1.Mean.Long   + Var.Y*rnorm(CUT2*N)
  Y.Z1.T0.Long <- Y.Z1.T0.Mean.Long   + Var.Y*rnorm(CUT2*N) 
  Y.Z1.T0.Long[rep(1:CUT2,N)<=k0] <- Y.Z1.T1.Long[rep(1:CUT2,N)<=k0]
  Y.Z0.T0.Long <- Y.Z0.T0.Mean.Long   + Var.Y*rnorm(CUT2*N)
  
  #############################
  # Potential Testing Horizon
  #############################
  
  H0 <- sapply(1:N,function(ii){H0.Ft(as.numeric(Y.Z0.T0.Long[k0+CUT2*(ii-1)]>=Yk0.threshold),
                                      X1[ii]+X2[ii])})
  H1 <- sapply(1:N,function(ii){H1.Ft(as.numeric(Y.Z1.T1.Long[k0+CUT2*(ii-1)]>=Yk0.threshold),
                                      X1[ii]+X2[ii])})
  H <- Z*H1 + (1-Z)*H0
  T <- as.numeric(H>CUT)
  TZ1 <- as.numeric(H1>CUT)
  Prop.H0 <- as.numeric(table(H0)/N)
  
  H0.Long <- rep(H0,each=CUT2)
  H1.Long <- rep(H1,each=CUT2)
  H.Long <- rep(H,each=CUT2)
  K.Long <- rep(1:CUT2,N)
  T.Long <- rep(T,each=CUT2)
  TZ1.Long <- rep(TZ1,each=CUT2)
  X1.Long <- rep(X1,each=CUT2)
  X2.Long <- rep(X2,each=CUT2)
  Y.Long <- Z.Long*T.Long*Y.Z1.T1.Long + Z.Long*(1-T.Long)*Y.Z1.T0.Long + (1-Z.Long)*Y.Z0.T0.Long
  Y.Z1.TZ.Long <- Y.Z1.T1.Long*TZ1.Long + Y.Z1.T0.Long*(1-TZ1.Long)
  Yk0.Long <- as.numeric( rep( Y.Long[k0+CUT2*((1:N)-1)] ,each=CUT2) >= Yk0.threshold )
  
  
  Data.Obs_all <- data.frame(ID=ID.Long,
                             Z=Z.Long,
                             H=H.Long,
                             H1=H1.Long,
                             H0=H0.Long,
                             K=K.Long,
                             T=T.Long,
                             X1=X1.Long,
                             X2=X2.Long,
                             
                             Y.Z1.T1=Y.Z1.T1.Mean.Long,
                             Y.Z1.TZ=Y.Z1.TZ.Mean.Long,
                             Y.Z1.T0=Y.Z1.T0.Mean.Long,
                             Y.Z0.T0=Y.Z0.T0.Mean.Long,
                             
                             Y.Z1.T1.Obs=Y.Z1.T1.Long,
                             Y.Z1.TZ.Obs=Y.Z1.TZ.Long,
                             Y.Z1.T0.Obs=Y.Z1.T0.Long,
                             Y.Z0.T0.Obs=Y.Z0.T0.Long,
                             
                             Y.Z1.T1.True=Y.Z1.T1.True,
                             Y.Z1.TZ.True=Y.Z1.TZ.True,
                             Y.Z1.T0.True=Y.Z1.T0.True,
                             Y.Z0.T0.True=Y.Z0.T0.True,
                             
                             Y=Y.Long,
                             Yk0=Yk0.Long)
  
  Data.Obs <- Data.Obs_all[Data.Obs_all$H>=Data.Obs_all$K,] # we cannot observe k larger than H. That is, cases with k <= H will be observed.
  
  return(list(Data.All=Data.Obs_all, Data.Obs=Data.Obs))
}

Y.Z0.B.Ft.Lin <- function(k,h,x,CUT){                                                        # expected potential outcome for whom used A=0 and H<=30 (i.e., non-users)
  CUT.T <- CUT
  t <- as.numeric(h-0.5>CUT)
  if(t==0){
    k <- min(k,CUT.T)
  }
  0.25*k+5
}

Y.Z1.B.Ft.Lin <- function(k,h,x,CUT){
  CUT.T <- CUT
  t <- as.numeric(h-0.5>CUT)
  if(t==0){
    k <- min(k,CUT.T)
    0.25*k+x*2+5
  } else {
    (0.25*k+x*2+5) + (k-k0)*0.1*as.numeric(k>k0)
  }
}

Y.Z0.B.Ft.Cub <- function(k,h,x,CUT){                                                        # expected potential outcome for whom used A=0 and H<=30 (i.e., non-users)
  CUT.T <- CUT
  t <- as.numeric(h-0.5>CUT)
  if(t==0){
    k <- min(k,CUT.T)
  }
  (-abs((k-1.25*CUT.T))^3+(1.25*CUT.T)^3)/((1.25*CUT.T)^3)*10 + x*2
}

Y.Z1.B.Ft.Cub <- function(k,h,x,CUT){
  CUT.T <- CUT
  t <- as.numeric(h-0.5>CUT)
  if(t==0){
    k <- min(k,CUT.T)
    (-abs((k-1.25*CUT.T))^3+(1.25*CUT.T)^3)/((1.25*CUT.T)^3)*(10) +2*x
  } else {
    ((-abs((k-1.25*CUT.T))^3+(1.25*CUT.T)^3)/((1.25*CUT.T)^3)*(10))*as.numeric(k<=1.25*CUT.T) + 
      10*as.numeric(k>1.25*CUT.T) + 2*x + 
      (-abs(((k-k0)-2.5*CUT.T))^3+(2.5*CUT.T)^3)/((2.5*CUT.T)^3)*(10)*as.numeric(k>k0)
  }
}






