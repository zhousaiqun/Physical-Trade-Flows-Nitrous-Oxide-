


load(str_c(pathout2,"/FABIO RegSec resolution.Rdata"))
load(str_c(pathout2,"/N2O satellite acc, GDP, Pop.Rdata"))
load(str_c(pathout3,"/N2O footprint of food consumption.Rdata"))

#PART 2: driving factors of total N2O emissions
#----------------------------------------------------
TarYear <- c(2020)
for(k in TarYear){
  YEAR<-c(2000,k)
  #YEAR<-c(k-1,k)#Chaining Decomposition
  pathout_sub <- paste(pathout4,"/",YEAR[1],"-",YEAR[2],sep="");dir.create(pathout_sub)#export detailed results
  print(paste(YEAR[1],"-",YEAR[2]))
  for(i in 1:length(YEAR)){
    load(str_c(pathout2,"/FABIO elements_",YEAR[i],".Rdata"))
    
    w<-as.vector(N2O_allyear[7,,which(years %in% YEAR[i])]/TotOutput)
    w[is.na(w)]<-0; w[w==Inf]<-0;w[w==-Inf]<-0
    
    FD_reg <- FD[,which(fd_codes$fd %in% "food")]
    
    FD_D <- array(0,dim = dim(FD_reg))
    
    for(m in 1:length(Regnam)){
      p=1+(m-1)*N;q=N+(m-1)*N
      FD_D[p:q,m]<-FD_reg[p:q,m]
    }
    FD_F <- FD_reg-FD_D
    FD_D_lev <- diag(colSums(FD_D));FD_F_lev <- diag(colSums(FD_F))
    FD_D_str <- t(t(FD_D)/colSums(FD_D));FD_F_str <- t(t(FD_F)/colSums(FD_F))
    FD_D_str[is.nan(FD_D_str)] <- 0;FD_F_str[is.nan(FD_F_str)] <- 0
    
    #separate trade and technology/patterns for intermediate and final demand
    fac <- rep(1:N,G)
    #A
    dd <- t(rowsum(as.matrix(A),fac,reorder = F))
    dim(dd) <- c(GN*N,1)
    Tech <- matrix(dd,nrow = GN,ncol = GN,byrow = T)
    Trd_int <- A/Tech
    Trd_int[is.nan(Trd_int)] <- 0
  
    if (i==1){
      FD_reg0<-FD_reg;FD_D_lev0 <- FD_D_lev;FD_D_str0 <- FD_D_str
      FD_F_lev0 <- FD_F_lev;FD_F_str0 <- FD_F_str;w0<-w;B0<-Leontief
      Tech0 <- Tech;Trd_int0 <- Trd_int
      rm(FD_reg,FD_D_lev,FD_D_str,FD_F_lev,FD_F_str,Leontief,TotOutput,Tech,Trd_int)
      gc()
    }else{
      FD_reg1<-FD_reg;FD_D_lev1 <- FD_D_lev;FD_D_str1 <- FD_D_str
      FD_F_lev1 <- FD_F_lev;FD_F_str1 <- FD_F_str;w1<-w;B1<-Leontief
      Tech1 <- Tech;Trd_int1 <- Trd_int
      rm(FD_reg,FD_D_lev,FD_D_str,FD_F_lev,FD_F_str,Leontief,TotOutput,Tech,Trd_int)
      gc()
    }
  }
  
  #decomposition
  #----------------------------------------
  tt <- Sys.time()
  Con.N2O <- N2Ofootprint[7,,]
  Chg_Con.N2O<-Con.N2O[,which(years %in% YEAR[2])]-Con.N2O[,which(years %in% YEAR[1])]
  SDA_Con.N2O<-array(0,dim = c(G,11),
                     dimnames = list(Regnam,
                                     c("Change in CB_N2O","E.Inten of Out_d","E.Inten of Out_f",
                                       "Cons.Str_d","Cons.Str_f","Cons.Lev_d","Cons.Lev_f",
                                       "Leontief_d_narrow_trade","Leontief_d_narrow_trade","Leontief_d_narrow_tech","Leontief_f_narrow_tech")))
  
  SDA_Con.N2O[,1] <- Chg_Con.N2O#Change in consumption-based emissions
  
  
  w_effect <- 0.5*(Matrix(diag(w1)-diag(w0)))%*%(B0%*%FD_reg0+B1%*%FD_reg1)
  
  Trd <- (0.5*B1%*%((Trd_int1-Trd_int0)*(Tech1+Tech0))%*%B0)
  B_trd_effect <- 0.5*(Matrix(diag(w0))%*%(Trd%*%FD_reg1)+diag(w1)%*%(Trd%*%FD_reg0))
  
  Tech <- (0.5*B1%*%((Trd_int0+Trd_int1)*(Tech1-Tech0))%*%B0)
  
  B_tech_effect <- 0.5*(Matrix(diag(w0))%*%(Tech%*%FD_reg1)+diag(w1)%*%(Tech%*%FD_reg0))
  
  w_effect_G<- array(0,dim=c(G,G));
  B_trd_effect_G<- array(0,dim=c(G,G))
  B_tech_effect_G<- array(0,dim=c(G,G))

  for(m in 1:length(Regnam)){
    p=1+(m-1)*N;q=N+(m-1)*N
    w_effect_G[m,]<-colSums(w_effect[p:q,])
    B_trd_effect_G[m,]<-colSums(B_trd_effect[p:q,])
    B_tech_effect_G[m,]<-colSums(B_tech_effect[p:q,])
  }
  
  SDA_Con.N2O[,2] <- diag(w_effect_G)[1:G]#domestic intensity effect
  SDA_Con.N2O[,3] <- colSums(w_effect_G)[1:G]-diag(w_effect_G)[1:G]#foreign intensity effect
  
  SDA_Con.N2O[,8] <- diag(B_trd_effect_G)[1:G]#domestic input trade effect
  SDA_Con.N2O[,9] <- colSums(B_trd_effect_G)[1:G]-diag(B_trd_effect_G)[1:G]#foreign input trade effect
  SDA_Con.N2O[,10] <- diag(B_tech_effect_G)[1:G]#domestic production technology effect
  SDA_Con.N2O[,11] <- colSums(B_tech_effect_G)[1:G]-diag(B_tech_effect_G)[1:G]#foreign production technology effect
  
  F_Str_effect_d <- 0.5*(Matrix(diag(w0))%*%(B0%*%((FD_D_str1-FD_D_str0)%*%FD_D_lev1))+
                           Matrix(diag(w1))%*%(B1%*%((FD_D_str1-FD_D_str0)%*%FD_D_lev0)))
  F_Str_effect_f <- 0.5*(Matrix(diag(w0))%*%(B0%*%((FD_F_str1-FD_F_str0)%*%FD_F_lev1))+
                           Matrix(diag(w1))%*%(B1%*%((FD_F_str1-FD_F_str0)%*%FD_F_lev0)))
  
  SDA_Con.N2O[,4] <- colSums(F_Str_effect_d)[1:G]
  SDA_Con.N2O[,5] <- colSums(F_Str_effect_f)[1:G]
  
  F_Lev_effect_d <- 0.5*(Matrix(diag(w1))%*%((B1%*%FD_D_str1)%*%(FD_D_lev1-FD_D_lev0))+
                                  Matrix(diag(w0))%*%((B0%*%FD_D_str0)%*%(FD_D_lev1-FD_D_lev0)))
  F_Lev_effect_f <- 0.5*(Matrix(diag(w1))%*%((B1%*%FD_F_str1)%*%(FD_F_lev1-FD_F_lev0))+
                                  Matrix(diag(w0))%*%((B0%*%FD_F_str0)%*%(FD_F_lev1-FD_F_lev0)))
  
  SDA_Con.N2O[,6] <- colSums(F_Lev_effect_d)[1:G]
  SDA_Con.N2O[,7] <- colSums(F_Lev_effect_f)[1:G]
  
  fnm = paste( pathout,"/",YEAR[1],"-",YEAR[2],"SDA_CBE_Effect.csv", sep="" );
  write.table(SDA_Con.N2O[,], file=fnm,sep="," )
  
  rownames(w_effect_G) <- Regnam;colnames(w_effect_G) <- Regnam
  
  dimnames(B_trd_effect_G) <- dimnames(w_effect_G)
  dimnames(B_tech_effect_G) <- dimnames(w_effect_G)
  
  fnm = paste( pathout_sub,"/",YEAR[1],"-",YEAR[2],"W_effect.csv", sep="" );
  write.table(w_effect_G[,], file=fnm,sep="," )
  fnm = paste( pathout_sub,"/",YEAR[1],"-",YEAR[2],"B_trd_effect.csv", sep="" );
  write.table(B_trd_effect_G[,], file=fnm,sep="," )
  fnm = paste( pathout_sub,"/",YEAR[1],"-",YEAR[2],"B_tech_effect.csv", sep="" );
  write.table(B_tech_effect_G[,], file=fnm,sep="," )
  #----------------------------------------
  
  rm(FD_reg0,FD_D_lev0,FD_D_str0,FD_F_lev0,FD_F_str0,w0,B0,
     FD_reg1,FD_D_lev1,FD_D_str1,FD_F_lev1,FD_F_str1,w1,B1,
     w_effect_G,w_effect,B_trd_effect_G,B_trd_effect,B_tech_effect_G,B_tech_effect,
     F_Lev_effect_d,F_Lev_effect_f,F_Str_effect_d,F_Str_effect_f)
  gc()
  print(paste(YEAR[1],"-",YEAR[2],"time cost: ",round(Sys.time( )-tt,3)) )
}
#----------------------------------------------------
rm(list = ls()[-which(ls() %in% c("path","pathout","pathout2","pathout3","pathout4",
                                  "pathdata","pathcode"))])
gc()
