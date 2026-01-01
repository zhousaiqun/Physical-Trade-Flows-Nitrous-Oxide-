
load(str_c(pathout2,"/FABIO RegSec resolution.Rdata"))
load(str_c(pathout2,"/N2O satellite acc, GDP, Pop.Rdata"))
Sec_group <- unique(io_codes$comm_group)

#focus on the N2O footprint of food consumption-----------
#country level
N2Ofootprint <- array(0,dim = c(length(SourceName),G,length(years)),
                      dimnames = list(SourceName,Regnam,years))#kt

Food_pattern <- array(0,dim = c(N,G,length(years)),
                      dimnames = list(Secnam,Regnam,years))

Food_pattern_Group <- array(0,dim = c(length(Sec_group),G,length(years)),
                            dimnames = list(Sec_group,Regnam,years))

N2O_FT_pattern <- Food_pattern
N2O_FT_pattern_Group <- Food_pattern_Group


Inten_pattern  <- Food_pattern
Inten_pattern_Group <- Food_pattern_Group

N2O_PerCap <- N2Ofootprint
N2O_PerCap_Group <- Food_pattern_Group
N2O_PerCap_Detail <- Food_pattern


#regional level
N2Ofootprint_agg <- array(0,dim = c(length(SourceName),length(REG18),length(years)),
                          dimnames = list(SourceName,REG18,years))#kt

Food_pattern_agg <- array(0,dim = c(length(Sec_group),length(REG18),length(years)),
                          dimnames = list(Sec_group,REG18,years))
N2O_FT_pattern_agg <- Food_pattern_agg
N2O_FT_pattern_agg_percap <- Food_pattern_agg
Inten_pattern_agg <-  Food_pattern_agg
N2O_PerCap_agg <- N2Ofootprint_agg

TradeShare <- rep(NA,length(years));names(TradeShare) <- years
TradeAmount <- TradeShare

Trade_reg_secgroup <- array(0,dim = c(3*length(Sec_group), G,length(years)),
                   dimnames = list(str_c(rep(c("Domestic","Import","Export"),each = length(Sec_group)),rep(Sec_group,3),sep = " : "),
                                   Regnam,years))#only calculate for total emission
Trade_reg_reg <- array(0,dim = c(G, G,length(years)),
                       dimnames = list(Regnam,Regnam,years))


fac <- io_codes$item
fac2 <- io_codes$comm_group
regfac <- regions$region

for (i in 1:length(years)) {
  t0 <- Sys.time()
  print(str_c("Begin calculation   ", years[i]))
  
  load(str_c(pathout2,"/FABIO elements_",years[i],".Rdata"))
  
  N_coef <- t(t(N2O_allyear[,,i])/TotOutput)
  N_coef[is.na(N_coef)] <- 0;N_coef[is.nan(N_coef)] <- 0
  N_coef[is.infinite(N_coef)] <- 0
  
  FoodDemand <- FD[,which(fd_codes$fd %in% "food")]
  Output_embodied <- Leontief%*%FoodDemand
  
  
  #N2O in trade, for total
  INT <- as.matrix(diag(N_coef[7,])%*%Output_embodied)
  INT2 <- INT
  for (k in 1:G) {
    m = N*(k-1)+1; n = N*k
    INT2[m:n,k] <- 0
  }
  
  Trade_reg_secgroup[1:length(Sec_group),,i] <- rowsum(INT-INT2,io_codes$comm_group,reorder = F)
  Trade_reg_secgroup[(length(Sec_group)+1):(2*length(Sec_group)),,i] <- rowsum(INT2,io_codes$comm_group,reorder = F)
  INT3 <- rowSums(INT2)
  dim(INT3) <- c(N,G)
  Trade_reg_secgroup[(2*length(Sec_group)+1):(3*length(Sec_group)),,i] <- rowsum(INT3,io_codes$comm_group[1:N],reorder = F)
  
  Trade_reg_reg[,,i] <- rowsum(INT,rep(Regnam,each = N),reorder = F)
  TradeShare[i] <- (sum(Trade_reg_reg[,,i],na.rm = T)-sum(diag(Trade_reg_reg[,,i]),na.rm = T))/sum(Trade_reg_reg[,,i],na.rm = T)
  TradeAmount[i] <- (sum(Trade_reg_reg[,,i],na.rm = T)-sum(diag(Trade_reg_reg[,,i]),na.rm = T))
  
   if (i == length(years)){
    IM_N2O_2020_reg <- rowsum(as.matrix((colSums(Trade_reg_reg[,,i],na.rm = T)-diag(Trade_reg_reg[,,i]))),
                              regfac,reorder = F,na.rm = T)/rowsum(as.matrix(colSums(Trade_reg_reg[,,i],na.rm = T)),
                                                                   regfac,reorder = F,na.rm = T)
  }
  
  for (k in 1:dim(N_coef)[1]) {
    N2Ofootprint[k,,i] <- as.vector(N_coef[k,]%*%Output_embodied)#national level
    N2O_PerCap[k,,i] <- N2Ofootprint[k,,i]/POP_allyear[,i]
    
    N2Ofootprint_agg[k,,i] <- rowsum(as.matrix(N2Ofootprint[k,,i]),regfac,reorder = F,na.rm = T)[-21]
    
    pop_year <- rowsum(as.matrix(POP_allyear[,i]),regfac,reorder = F,na.rm = T)[-21]
    N2O_PerCap_agg[k,,i] <- N2Ofootprint_agg[k,,i]/pop_year
    
    if (k == 7){
      N_multi <- as.vector(N_coef[k,]%*%Leontief)
      
      INT <- array(0,dim = c(GN,G))
      for (r in 1:G) {
        INT[,r] <- N_multi*FoodDemand[,r]
      }
      
      N2O_FT_pattern[,,i] <- rowsum(INT,fac,reorder = F,na.rm = T)
      N2O_FT_pattern_Group[,,i] <- rowsum(INT,fac2,reorder = F,na.rm = T)
      N2O_FT_pattern_agg[,,i] <- t(rowsum(t(N2O_FT_pattern_Group[,,i]),regfac,reorder = F,na.rm = T))[,-21]
      N2O_FT_pattern_agg_percap[,,i] <- t(t(N2O_FT_pattern_agg[,,i])/pop_year)
      
      N2O_PerCap_Detail[,,i] <- N2O_FT_pattern[,,i]/POP_allyear[,i]
      N2O_PerCap_Group[,,i] <- N2O_FT_pattern_Group[,,i]/POP_allyear[,i]
      
    }
  }
  
  #Food consumption pattern
  Food_pattern[,,i] <- rowsum(as.matrix(FoodDemand),fac,reorder = F,na.rm = T)
  Food_pattern_Group[,,i] <- rowsum(as.matrix(FoodDemand),fac2,reorder = F,na.rm = T)
  Food_pattern_agg[,,i] <- t(rowsum(t(Food_pattern_Group[,,i]),regfac,reorder = F,na.rm = T))[,-21]
  
  #Intensity pattern
  Inten_pattern[,,i] <- N2O_FT_pattern[,,i]/Food_pattern[,,i]
  Inten_pattern_Group[,,i] <- N2O_FT_pattern_Group[,,i]/Food_pattern_Group[,,i]
  Inten_pattern_agg[,,i] <- t(rowsum(t(Inten_pattern_Group[,,i]),regfac,reorder = F,na.rm = T))[,-21]
  
  #country average Intensity pattern
  
  In <- matrix(N_coef[7,], nrow = 1) %*% Leontief %*% diag(rowSums(FoodDemand))
  N_av_intensity <- rowsum(as.numeric(In), rep(Regnam, each = N), reorder = FALSE) / 
    rowsum(as.numeric(rowSums(FoodDemand)), rep(Regnam, each = N), reorder = FALSE)
  In_total <- sum(In) 
  Food_total <- sum(FoodDemand)
  N_global_intensity <- sum(In) / sum(FoodDemand)
  write.csv(N_av_intensity, str_c(pathout3, "/ave_N2O_Intensity_perunit_of_FoodConsumption_", years[i], ".csv"))
  write.csv(data.frame(Global_N2O_Intensity = N_global_intensity),
            str_c(pathout3, "/global_N2O_Intensity_", years[i], ".csv"))
  

  # 计算每个国家的N₂O总嵌入量
  In_total_by_country <- rowsum(as.numeric(In), rep(Regnam, each = N), reorder = FALSE)
  
  # 计算每个国家的总食物消费量
  Food_total_by_country <- rowsum(as.numeric(rowSums(FoodDemand)), rep(Regnam, each = N), reorder = FALSE)
  
  # 合并为一个数据框（不计算强度）
  Both_df <- data.frame(
    Country = rownames(Food_total_by_country),
    N2O_Embedded_kg = as.numeric(In_total_by_country),
    FoodConsumption_kg = as.numeric(Food_total_by_country)
  )
  
  # 写入CSV（不含N2O强度）
  write.csv(Both_df, str_c(pathout3, "/N2O_and_Food_total_", years[i], ".csv"), row.names = FALSE)
  
  #N2O in trade_sector level (from demand perspective), add on July 15, 2025======
 if(years[i] %in% c(2020)){
   Trade_reg_secgroup_con <- array(0,dim = c(3*length(Sec_group), G),
                                   dimnames = list(str_c(rep(c("Domestic","Import","Export"),each = length(Sec_group)),rep(Sec_group,3),sep = " : "),
                                                   Regnam))
   
   Embodied_inten <- diag(N_coef[7,])%*%Leontief
   fac3 <- rep(regions$area,each = N)
   Embodied_inten_reg <- rowsum(as.matrix(Embodied_inten),fac3,reorder = F)
   
   INT4 <- array(0, dim = c(G*N,G))
   for (r in 1:G) {
     print(str_c("Begin  ",Regnam[r]))
     test <- rowsum(t(Embodied_inten_reg%*%diag(FoodDemand[,r])),fac,reorder = F)
     dim(test) <- c(G*N,1)
     INT4[,r] <- test
   }
   
   INT5 <- INT4
   for (k in 1:G) {
     m = N*(k-1)+1; n = N*k
     INT5[m:n,k] <- 0
   }
   
   Trade_reg_secgroup_con[1:length(Sec_group),] <- rowsum(INT4-INT5,io_codes$comm_group,reorder = F)
   Trade_reg_secgroup_con[(length(Sec_group)+1):(2*length(Sec_group)),] <- rowsum(INT5,io_codes$comm_group,reorder = F)
   INT6 <- rowSums(INT5)
   dim(INT6) <- c(N,G)
   Trade_reg_secgroup_con[(2*length(Sec_group)+1):(3*length(Sec_group)),] <- rowsum(INT6,io_codes$comm_group[1:N],reorder = F)
   
   write.csv(Trade_reg_secgroup_con,
             file = str_c(pathout3,"/N2O trade_domestic, import and export_cons perspective",years[i],".csv"))
   
 }
  #======
  
  #save as csv file
  write.csv(N2Ofootprint[,,i],file = str_c(pathout3,"/N2O footprint by source_national",years[i],".csv"))
  write.csv(N2Ofootprint_agg[,,i],file = str_c(pathout3,"/N2O footprint by source_regional",years[i],".csv"))

  write.csv(Food_pattern[,,i],file = str_c(pathout3,"/Food consumption pattern_national_123 prods",years[i],".csv"))
  write.csv(N2O_FT_pattern[,,i],file = str_c(pathout3,"/N2O footprint pattern_national_123 prods",years[i],".csv"))
  write.csv(Inten_pattern[,,i],file = str_c(pathout3,"/N2O intensity pattern_national_123 prods",years[i],".csv"))

  write.csv(Food_pattern_Group[,,i],file = str_c(pathout3,"/Food consumption pattern_national_23 groups",years[i],".csv"))
  write.csv(N2O_FT_pattern_Group[,,i],file = str_c(pathout3,"/N2O footprint pattern_national_23 groups",years[i],".csv"))
  write.csv(Inten_pattern_Group[,,i],file = str_c(pathout3,"/N2O intensity pattern_national_23 groups",years[i],".csv"))

  write.csv(Food_pattern_agg[,,i],file = str_c(pathout3,"/Food consumption pattern_regional_23 groups",years[i],".csv"))
  write.csv(N2O_FT_pattern_agg[,,i],file = str_c(pathout3,"/N2O footprint pattern_regional_23 groups",years[i],".csv"))
  write.csv(N2O_FT_pattern_agg_percap[,,i],file = str_c(pathout3,"/N2O footprint pattern_PerCapita_regional_23 groups",years[i],".csv"))
  write.csv(Inten_pattern_agg[,,i],file = str_c(pathout3,"/N2O intensity pattern_regional_23 groups",years[i],".csv"))
  write.csv(N2O_PerCap_Detail[,,i],file = str_c(pathout3,"/N2O footprint pattern_PerCapita_national_123 prods",years[i],".csv"))
  write.csv(N2O_PerCap_Group[,,i],file = str_c(pathout3,"/N2O footprint pattern_PerCapita_national_groups",years[i],".csv"))
  write.csv(N2O_PerCap[,,i],file = str_c(pathout3,"/N2O footprint pattern_PerCapita_national_source",years[i],".csv"))
  
  write.csv(Trade_reg_reg[,,i],file = str_c(pathout3,"/N2O flow between countries",years[i],".csv"))
  write.csv(Trade_reg_secgroup[,,i],file = str_c(pathout3,"/N2O trade_domestic, import and export",years[i],".csv"))
  
  print(str_c("Time cost,", years[i],"------",round(Sys.time()-t0,2)))
}

write.csv(TradeShare,file = str_c(pathout3,"/TradeShare_global total_time series.csv"))
write.csv(TradeAmount,file = str_c(pathout3,"/TradeAmount (kt)_global total_time series.csv"))
write.csv(IM_N2O_2020_reg,file = str_c(pathout3,"/Import share by region_2020.csv"))

save(N2Ofootprint,Food_pattern,Food_pattern_Group,
     regfac,REG18,regions, 
     N2O_PerCap_agg,N2O_PerCap,
     N2O_FT_pattern,N2O_FT_pattern_Group,Inten_pattern,Inten_pattern_Group,
     N2Ofootprint_agg,Food_pattern_agg,N2O_FT_pattern_agg,Inten_pattern_agg,
     TradeShare,IM_N2O_2020_reg,
     file = str_c(pathout3,"/N2O footprint of food consumption.Rdata"))

rm(list = ls()[-which(ls() %in% c("path","pathout","pathout2","pathout3","pathout4",
                                  "pathdata","pathcode"))])
gc()
