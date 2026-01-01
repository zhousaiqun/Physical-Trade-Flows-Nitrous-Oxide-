

library(stringr)
library(data.table)
library(Matrix)


setwd("D:\\FABIO N2O_250525")
path <- getwd()
pathdata <- str_c(path,"/Data")
pathcode <- str_c(path,"/Code")
pathout <- str_c(path,"/Analysis_250525");dir.create(pathout)
pathout2 <- str_c(pathout,"/FABIO_Rdata");dir.create(pathout2)
pathout3 <- str_c(pathout,"/Final results");dir.create(pathout3)
pathout4 <- str_c(pathout,"/SDA results");dir.create(pathout4)


source(str_c(pathcode, "/Module 1_convert FABIO ras data to R data.R"))#Module 1: convert FABIO ras data to R data
source(str_c(pathcode, "/Module 2_clean the N2O satellite account.R"))#Module 2: clean the NO2 satellite account
source(str_c(pathcode, "/Module 3_Calculate N2O footprint_v2.R"))#Module 3: Calculate NO2 footprint at country and sectoral level
source(str_c(pathcode, "/Module 4_Driving factor analysis_two polar.R"))#Module 4: Driving factor analysis

