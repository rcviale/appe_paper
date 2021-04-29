# Load packages
library(rugarch)
library(quantmod)
library(forecast)

# File path (Rodrigo)
path <- 'C:\\Users\\rodri\\OneDrive\\Documents\\Academics\\Univerzita Karlova\\2nd semester\\Applied Econometrics\\^GSPC.csv'

# File path (Dominik)
#path <- '(paste your location for the data here)\\^GSPC.csv'

# Load data
raw_data <- read.csv(path)

# Compiute log returns
rets <- diff(log(raw_data[, 2]))

# ACF and PACF
acf(rets)
pacf(rets)
# Strong autocorrelations in lags 1 and 2

# Testing auto.arima function
autofit <- auto.arima(rets)
autofit
# Yields ARIMA(4, 0, 0), quite a lot of coefficients

# Estimating an AR(2) and checking BIC for both fits
ar2 <- arima(rets, order = c(2, 0, 0))
BIC(ar2)
BIC(autofit)
# BIC for AR(2) is better