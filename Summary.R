Result <- list()
Result[[1]] <- Result[[2]] <- Result[[3]] <- Result[[4]] <- Result[[5]] <- list()
Result.Cubic <- list()
Result.Cubic[[1]] <- Result.Cubic[[2]] <- Result.Cubic[[3]] <- Result.Cubic[[4]] <- Result.Cubic[[5]] <- list()
Folder <- "Result_F"
N.grid <- c(100,200,500,1000)

NT <- 4

for(nn in 1:NT){
  Result[[1]][[nn]] <- read.csv(sprintf("%s/Merge_NDE_Type1_N_%0.4d.csv",Folder,N.grid[nn]))
  Result[[2]][[nn]] <- read.csv(sprintf("%s/Merge_NIE_Type1_N_%0.4d.csv",Folder,N.grid[nn]))
  Result[[3]][[nn]] <- read.csv(sprintf("%s/Merge_TE_Type1_N_%0.4d.csv",Folder,N.grid[nn]))
  Result[[4]][[nn]] <- read.csv(sprintf("%s/Merge_Cover_Type1_N_%0.4d.csv",Folder,N.grid[nn]))
  
  Result.Cubic[[1]][[nn]] <- read.csv(sprintf("%s/Merge_NDE_Type2_N_%0.4d.csv",Folder,N.grid[nn]))
  Result.Cubic[[2]][[nn]] <- read.csv(sprintf("%s/Merge_NIE_Type2_N_%0.4d.csv",Folder,N.grid[nn]))
  Result.Cubic[[3]][[nn]] <- read.csv(sprintf("%s/Merge_TE_Type2_N_%0.4d.csv",Folder,N.grid[nn]))
  Result.Cubic[[4]][[nn]] <- read.csv(sprintf("%s/Merge_Cover_Type2_N_%0.4d.csv",Folder,N.grid[nn]))
  
  Result[[1]][[nn]] <- Result[[1]][[nn]][which(!is.na(apply(Result[[1]][[nn]],1,sum))),]
  Result[[2]][[nn]] <- Result[[2]][[nn]][which(!is.na(apply(Result[[2]][[nn]],1,sum))),]
  Result[[3]][[nn]] <- Result[[3]][[nn]][which(!is.na(apply(Result[[3]][[nn]],1,sum))),]
  Result[[4]][[nn]] <- Result[[4]][[nn]][which(!is.na(apply(Result[[4]][[nn]],1,sum))),]
  
  Result.Cubic[[1]][[nn]] <- Result.Cubic[[1]][[nn]][which(!is.na(apply(Result.Cubic[[1]][[nn]],1,sum))),]
  Result.Cubic[[2]][[nn]] <- Result.Cubic[[2]][[nn]][which(!is.na(apply(Result.Cubic[[2]][[nn]],1,sum))),]
  Result.Cubic[[3]][[nn]] <- Result.Cubic[[3]][[nn]][which(!is.na(apply(Result.Cubic[[3]][[nn]],1,sum))),]
  Result.Cubic[[4]][[nn]] <- Result.Cubic[[4]][[nn]][which(!is.na(apply(Result.Cubic[[4]][[nn]],1,sum))),]
  
}
CUT1 <- 20
CUT2 <- 40
k0 <- 5

BIAS <- function(MAT){
  BIAS.GAM <- MAT[,3*CUT2+(1+k0):CUT2]
  BIAS.Lin <- MAT[,1*CUT2+(1+k0):CUT2]
  BIAS.Cubic <- MAT[,2*CUT2+(1+k0):CUT2]
  BIAS.Oracle <- MAT[,0*CUT2+(1+k0):CUT2]
  c( mean(apply(BIAS.GAM,2,mean)),
     mean(apply(BIAS.Lin,2,mean)),
     mean(apply(BIAS.Cubic,2,mean)),
     mean(apply(BIAS.Oracle,2,mean)))
}

ESE <- function(MAT){
  BIAS.GAM <- MAT[,3*CUT2+(1+k0):CUT2]
  BIAS.Lin <- MAT[,1*CUT2+(1+k0):CUT2]
  BIAS.Cubic <- MAT[,2*CUT2+(1+k0):CUT2]
  BIAS.Oracle <- MAT[,0*CUT2+(1+k0):CUT2]
  c( mean(apply(BIAS.GAM,2,sd)),
     mean(apply(BIAS.Lin,2,sd)),
     mean(apply(BIAS.Cubic,2,sd)),
     mean(apply(BIAS.Oracle,2,sd)))
}

MSE <- function(MAT){
  BIAS.GAM <- MAT[,3*CUT2+(1+k0):CUT2]
  BIAS.Lin <- MAT[,1*CUT2+(1+k0):CUT2]
  BIAS.Cubic <- MAT[,2*CUT2+(1+k0):CUT2]
  BIAS.Oracle <- MAT[,0*CUT2+(1+k0):CUT2]
  c( mean(apply(BIAS.GAM^2,2,mean)),
     mean(apply(BIAS.Lin^2,2,mean)),
     mean(apply(BIAS.Cubic^2,2,mean)),
     mean(apply(BIAS.Oracle^2,2,mean)))
}

COVER <- function(MAT,
                  cover="Cover",
                  estimand="NDE"){
  pos <- which(substr(colnames(MAT),1,nchar(cover)+nchar(estimand)+1)==sprintf("%s_%s",cover,estimand))
  COVER.GAM <- MAT[,pos]
  c( apply(COVER.GAM,2,mean))
} 

Summary <- function(MAT){
  cbind(BIAS(MAT),
        ESE(MAT),
        MSE(MAT))
}

SUM <- list()
for(tt in 1:NT){
  SUM[[tt]] <- cbind( Summary(Result[[1]][[tt]]), 
                      Summary(Result[[2]][[tt]]),
                      Summary(Result[[3]][[tt]])) 
}

SSS <- rbind(SUM[[1]],SUM[[2]])
for(tt in 3:NT){
  SSS <- rbind(SSS,SUM[[tt]])
}

SUM.All <- cbind(rep(N.grid,each=4),
                 SSS)

colnames(SUM.All) <- c("N",
                       sprintf("NDE_%s",c("Bias","ESE","MSE")),
                       sprintf("NIE_%s",c("Bias","ESE","MSE")),
                       sprintf("TE_%s",c("Bias","ESE","MSE")))



SUM.Cubic <- list()
for(tt in 1:NT){
  SUM.Cubic[[tt]] <- cbind( Summary(Result.Cubic[[1]][[tt]]), 
                            Summary(Result.Cubic[[2]][[tt]]),
                            Summary(Result.Cubic[[3]][[tt]])) 
}

SSS.Cubic <- rbind(SUM.Cubic[[1]],SUM.Cubic[[2]])
for(tt in 3:NT){
  SSS.Cubic <- rbind(SSS.Cubic,SUM.Cubic[[tt]])
}

SUM.All.Cubic <- cbind(rep(N.grid,each=4),
                       SSS.Cubic)

colnames(SUM.All.Cubic) <- c("N",
                             sprintf("NDE_%s",c("Bias","ESE","MSE")),
                             sprintf("NIE_%s",c("Bias","ESE","MSE")),
                             sprintf("TE_%s",c("Bias","ESE","MSE")))


write.csv(SUM.All[-(1:4)*4,c(1,2,4,5,7,8,10)],
          "Summary_Lin.csv",
          row.names=F)

write.csv(SUM.All.Cubic[-(1:4)*4,c(1,2,4,5,7,8,10)],
          "Summary_Cubic.csv",
          row.names=F)

SUM.All[-(1:4)*4,c(1,2,4,5,7,8,10)]

SUM.All.Cubic[-(1:4)*4,c(1,2,4,5,7,8,10)]



CL <- colnames(Result[[1]][[tt]])

MAR <- c(3,3,2,1)
Valid <- (6):CUT2
Valid2 <- which(as.numeric(sapply(CL,
                 function(v){substr(v,nchar(v)-3,nchar(v))}))>5)

for(cover in c("Cover")){
  
  png(sprintf("%s_%s_Coverage.png",Folder,cover),height=4,width=8,unit="in",res=500)
  
  layout(matrix(c(19,20,13,14,15,16,18,
                  12,10,1,2,3,4,5,
                  12,11,6,7,8,9,5,
                  21,22,17,17,17,17,23),4,7,byrow=T),
         widths=c(1,1,5,5,5,5,2),
         heights=c(1,10,10,1))
  
  par(mar=MAR)
  
  for(tt in 1:NT){
    plot(Valid-0.2,
         COVER(Result[[4]][[tt]][,Valid2],cover=cover,estimand="NDE"),
         xlim=c(0,CUT2+1),
         ylim=c(0.9,1),
         pch=19,col=1,cex=0.5,
         yaxt="n",
         xlab="",
         ylab="")
    points(Valid,
           COVER(Result[[4]][[tt]][,Valid2],cover=cover,estimand="NIE"),
           pch=3,col=2,lwd=1)
    points(Valid+0.2,
           COVER(Result[[4]][[tt]][,Valid2],cover=cover,estimand="TE"),
           pch=4,col=4,lwd=1)
    abline(h=0.95,lty=2)
    axis(2,at=c(0.9,0.95,1))
    # title(main=sprintf("N=%s",N.grid[tt]))
  }
  
  par(mar=c(1,0,1,0)*MAR)
  
  plot.new()
  points(0.1,0.9,pch=19,col=1,cex=1.5)
  text(0.2,0.9,"NDE",pos=4)
  points(0.1,0.7,pch=3,col=2,cex=2,lwd=2)
  text(0.2,0.7,"NIE",pos=4)
  points(0.1,0.5,pch=4,col=4,cex=2,lwd=2)
  text(0.2,0.5,"TE",pos=4)
  
  par(mar=MAR)
  
  for(tt in 1:NT){
    plot(Valid-0.2,
         COVER(Result.Cubic[[4]][[tt]][,Valid2],cover=cover,estimand="NDE"),
         xlim=c(0,CUT2+1),
         ylim=c(0.9,1),
         pch=19,col=1,cex=0.5,
         yaxt="n",
         xlab="",
         ylab="")
    points(Valid,
           COVER(Result.Cubic[[4]][[tt]][,Valid2],cover=cover,estimand="NIE"),
           pch=3,col=2,lwd=1)
    points(Valid+0.2,
           COVER(Result.Cubic[[4]][[tt]][,Valid2],cover=cover,estimand="TE"),
           pch=4,col=4,lwd=1)
    abline(h=0.95,lty=2)
    axis(2,at=c(0.9,0.95,1))
    # title(main=sprintf("N=%s",N.grid[tt]))
  }
  
  par(mar=c(1,0,1,0)*MAR)
  plot.new()
  text(0.5,0.5,"Linear Outcome Model",srt=90,cex=1)
  plot.new()
  text(0.5,0.5,"Cubic Outcome Model",srt=90,cex=1)
  plot.new()
  text(0.5,0.5,"Coverage",srt=90,cex=1.5)
  
  par(mar=c(0,1,0,1)*MAR)
  for(nn in 1:4){
    plot.new()
    text(0.5,0.5,sprintf("N=%s",
                         N.grid[nn]),
         cex=1.25) 
  }
  
  
  par(mar=c(0,1,0,1)*MAR)
  plot.new()
  text(0.5,0.5,"Time",cex=1.5)
  plot.new()
  
  dev.off()
  
}
