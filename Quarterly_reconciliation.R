## Load libraries
library(hts)

## Set working directory to file location and load utils
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("matriz_y.RData")
load("hierarchy.RData")

## Transform monthly time series to the desired frequency
periodicity= 3
vetor_y    = sapply(seq(1,204,by=3),function(i) rowSums(vetor_y[,i:(i+periodicity-1)]))

## Forecast settings
# Load predicted values
y_hat  = read.csv('base_forecast/base_forecast_quarter.csv',header=FALSE,sep=';')
# Set the number of steps ahead to forecast
h      = 1
# Set the sample size for model training
tam_treino = 9*4
# Set the sample size for model testing
tam_teste  = dim(vetor_y)[2] - tam_treino

## Run reconciliation
recon_h1_ols  = NULL
recon_h1_atha = NULL
recon_h1_wls  = NULL
recon_h1_mint = NULL
for(i in 1:tam_teste){
  m_res    = read.csv(paste('residuals/Quarter/residuo',(tam_treino+i-1),'.csv',sep=''),header=FALSE,sep=';')
  grupos   = gts(y = t(vetor_y[87:271,1:(tam_treino+i-1)]),groups = m_groups)
  matriz_s = smatrix(grupos)
  
  # WLS_s covariance matrix
  aux_atha   = c(matriz_s%*%matrix(1,185,1))
  recon_atha = combinef(fcasts = matrix(y_hat[,i],1,byrow=TRUE), groups = grupos$groups, 
                        weights = aux_atha, nonnegative = TRUE, keep = 'all', algorithms = "lu")
  
  # WLS_upsilon covariance matrix
  aux_wls    = diag(as.matrix(m_res)%*%t(as.matrix(m_res)))/(tam_treino+i-1)
  recon_wls  = combinef(fcasts = matrix(y_hat[,i],1,byrow=TRUE), groups = grupos$groups, 
                        weights = aux_wls, nonnegative = TRUE, keep = 'all', algorithms = "lu")
  
  # MinT covariance matrix
  recon_mint = MinT(fcasts = matrix(y_hat[,i],1,byrow=TRUE), groups = grupos$groups, 
                    residual = t(m_res), covariance = "shr", nonnegative = TRUE, keep = 'all', 
                    algorithms = "lu")
  
  # OLS covariance matrix
  recon_ols = combinef(fcasts = matrix(y_hat[,i],1,byrow=TRUE), groups = grupos$groups, 
                        weights = rep(1,271), nonnegative = TRUE, keep = 'all', algorithms = "lu")
  
  recon_h1_atha = cbind(recon_h1_atha,t(recon_atha))
  recon_h1_wls  = cbind(recon_h1_wls,t(recon_wls))
  recon_h1_mint = cbind(recon_h1_mint,t(recon_mint))
  recon_h1_ols  = cbind(recon_h1_ols,t(recon_ols))
  
  print(i)
}