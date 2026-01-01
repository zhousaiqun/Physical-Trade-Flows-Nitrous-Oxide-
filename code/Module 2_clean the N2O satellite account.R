

#Need further take the emission from energy into account



load(str_c(pathout2,"/FABIO RegSec resolution.Rdata"))

#Read the N2O data-----
Test <- readxl::read_xlsx(str_c(pathdata,"/N2O_Emission_data/FAO_crop&livestock_emission_2000.xlsx"))
SourceName <- colnames(Test)[9:15]

N2O_allyear <- array(0,dim = c(length(SourceName),GN,length(years)),
                     dimnames = list(SourceName,Regsecnam,years))

for (i in 1:length(years)) {
  N2O_allyear[,,i] <- t(apply(readxl::read_xlsx(str_c(pathdata,"/N2O_Emission_data/FAO_crop&livestock_emission_",
                                                      years[i],".xlsx"))[,9:15], 2, as.matrix))
  
}
N2O_allyear[is.na(N2O_allyear)] <- 0



#GDP data---
RegMap_WB <- read.csv(str_c(pathdata,"/Region map_FABIO and WB.csv"))

GDP_allyear <- array(0,dim = c(G,length(years)),
                     dimnames = list(Regnam,years))#2015 constant price $

GDP_rawdata <- read.csv(str_c(pathdata,"/GDP_WORLD_BANK.csv"))
GDP_filter <- GDP_rawdata[match(RegMap_WB$ISO_WB[1:184],GDP_rawdata$Short.name),]
GDP_match <- GDP_filter[match(regions$iso3c,RegMap_WB$ISO3[1:184]),]

GDP_allyear[,] <- apply(GDP_match[,which(colnames(GDP_match) %in% "X2000"):which(colnames(GDP_match) %in% "X2020")], 2, as.numeric)




#Population data-----
POP_allyear <- array(0,dim = c(G,length(years)),
                     dimnames = list(Regnam,years))#headcount

POP_rawdata <- read.csv(str_c(pathdata,"/Population_WORLD_BANK.csv"))
POP_filter <- POP_rawdata[match(RegMap_WB$ISO_WB[1:184],POP_rawdata$Short.name),]
POP_match <- POP_filter[match(regions$iso3c,RegMap_WB$ISO3[1:184]),]

POP_allyear[,] <- apply(POP_match[,which(colnames(POP_match) %in% "X2000"):which(colnames(POP_match) %in% "X2020")], 2, as.numeric)


#Save data-----
save(N2O_allyear,POP_allyear,GDP_allyear, SourceName,
     file = str_c(pathout2,"/N2O satellite acc, GDP, Pop.Rdata"))

rm(list = ls()[-which(ls() %in% c("path","pathout","pathout2","pathout3",
                                  "pathdata","pathcode"))])
gc()
