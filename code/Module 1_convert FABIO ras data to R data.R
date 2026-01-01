

#The function for calculating Leontief inverse is obtained fr https://github.com/fineprint-global/fabio/blob/master/R/11_leontief_inverse.R



# Leontief inverse function ---
prep_solve <- function(year, Z, X,
                       adj_X = FALSE, adj_A = TRUE, adj_diag = FALSE) {
  
  if(adj_X) {X <- X + 1e-10}
  
  A <- Matrix(0, nrow(Z), ncol(Z))
  idx <- X != 0
  A[, idx] <- t(t(Z[, idx]) / X[idx])
  #A <- Z
  #A@x <- A@x / rep.int(X, diff(A@p))
  
  if(adj_A) {A[A < 0] <- 0}
  if(adj_diag) {diag(A)[diag(A) == 1] <- 1 - 1e-10}
  
  L <- .sparseDiagonal(nrow(A)) - A
  
  lu(L) # Computes LU decomposition and stores it in L
  
  #tryCatch({
  L_inv <- solve(L, tol = .Machine[["double.eps"]], sparse = TRUE)
  #}, error=function(e){cat("ERROR in ", year, "\n")})
  
  #L_inv[L_inv<0] <- 0
  
  return(L_inv)
}


# Loop for years to save the basic MRIO element------------
years <- 2000:2020
years_singular_losses <- c(2006, 2007, 2010, 2011, 2013, 2018, 2019) #  c(2013,2019) #c(1990,2010,2019) #c(1994,2002,2009)


# Leontief inverse for losses version of fabio
X <- readRDS(str_c(pathdata,"/fabio_v1.2_losses/X.rds"))
Y <- readRDS(str_c(pathdata,"/fabio_v1.2_losses/Y.rds"))
Z_v <- readRDS(str_c(pathdata,"/fabio_v1.2_losses/Z_value.rds"))
E_data <- readRDS(str_c(pathdata,"/fabio_v1.2_losses/E.rds"))

for(year in years){
  t0 <- Sys.time()
  print(str_c("Begin calculation   ", year))
  
  
  Z <- Z_v[[as.character(year)]]
  X_filter <-  X[, as.character(year)]
  A <- Matrix(0, nrow(Z), ncol(Z))
  
  idx <- X_filter != 0
  A[, idx] <- t(t(Z[, idx]) / X_filter[idx])
  
  adjust_losses <- ifelse(year %in% years_singular_losses, TRUE, FALSE)
  
  Leontief <- prep_solve(year = year, Z = Z_v[[as.character(year)]],
                  X = X[, as.character(year)], adj_diag = adjust_losses)
  Leontief[Leontief<0] <- 0
  
  FD <- Y[[as.character(year)]]
  TotOutput <- X[, as.character(year)]
  InterTrade <- Z_v[[as.character(year)]]
  
  Envior <- E_data[[as.character(year)]]
  
  save(FD,TotOutput,InterTrade,Leontief,Envior,A,
       file = str_c(pathout2,"/FABIO elements_",year,".Rdata"))
  
  print(str_c("Time cost,", year,"------",round(Sys.time()-t0,2)))
}


#Load and save the regnam, secnam, envir satellite acccount-----
fd_codes <- read.csv(str_c(pathdata,"/fabio_v1.2_losses/fd_codes.csv"))
io_codes <- read.csv(str_c(pathdata,"/fabio_v1.2_losses/io_codes.csv"))

items <- read.csv(str_c(pathdata,"/fabio_v1.2_losses/items.csv"))
regions <- read.csv(str_c(pathdata,"/fabio_v1.2_losses/regions_v2.csv"))

Regnam <- as.character(regions$area)
G <- length(Regnam)
Secnam <- as.character(items$item)
N <- length(Secnam)
Regsecnam <- str_c(rep(Regnam, each = N),rep(Secnam, G))
GN <- length(Regsecnam)

io_codes_new <- read.csv(str_c(pathdata,"/group_codes.csv"))
io_codes$comm_group <- rep(io_codes_new$comm_group_V2,G)

REG18 <- unique(regions$region)[-21]

save(fd_codes,io_codes,items,regions,
     Regnam,Secnam,Regsecnam,G,N,GN,years,REG18,
     file = str_c(pathout2,"/FABIO RegSec resolution.Rdata"))

rm(list = ls()[-which(ls() %in% c("path","pathout","pathout2","pathout3","pathout4",
                                  "pathdata","pathcode"))])
gc()
