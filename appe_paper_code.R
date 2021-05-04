# Load packages
library(rugarch)
library(quantmod)
library(forecast)
library(tseries)
library(aTSA)

# File path (Rodrigo)
path <- 'C:\\Users\\rodri\\OneDrive\\Documents\\Academics\\Univerzita Karlova\\2nd semester\\Applied Econometrics\\^GSPC.csv'

# File path (Dominik)
#path <- '(paste your location for the data here)\\^GSPC.csv'

# Load data
raw_data <- read.csv(path)

# Plot for price level
plot.ts(raw_data[, 2], main = 'GSPC Prices', ylab = NA)

# Test for stationarity
adf.test(raw_data[, 2])
# Rejects alternative of stationarity

# Compute log returns
rets <- diff(log(raw_data[, 2]))

# Plot for log returns
plot.ts(rets, main = 'Log Returns', ylab = NA)

# Testing for stationarity
adf.test(rets)
# Rejects null of unit root

# ACF and PACF
acf(rets)
pacf(rets)
# Strong autocorrelations in lags 1 and 2 suggests we should model mean (ARMA)

# Testing auto.arima function suggestion
autofit <- auto.arima(rets)
autofit
#   Yields ARIMA(4, 0, 0), quite a lot of coefficients, considering we will have
# GARCH coefficients after that

# Estimating other ARMA and checking  BIC for both fits
crits <- matrix(ncol = 2, nrow = 25)
colnames(crits) <- c('Model', 'BIC')
for (i in 1 : 5){
  for (j in 1 : 5){
    model <- arima(rets, order = c((i - 1), 0, (j - 1)))
    crits[5 * (i - 1) + j, 1] <- paste0('ARMA(', (i - 1), ', ', (j - 1), ')')
    crits[5 * (i - 1) + j, 2] <- BIC(model)
  }
}
View(crits)

# To avoid overparametrization, we will use AR(2) based on BIC and parsimonity
arma43 <- arima(rets, order = c(4, 0, 3))

rm(crits, i, j, model, path)

# Test for ARCH effects
summary(arch.test(arma43))
# Rejects null of homoskedasticity

rm(arma43)

# Trying GARCH(1, 1) with the 5 best BIC ARMA models
sg11_43 <- ugarchspec(mean.model = list(armaOrder = c(4, 3)),
                      variance.model = list(garchOrder = c(1, 1), model = "sGARCH"),
                      distribution.model = 'norm')
sg11_44 <- ugarchspec(mean.model = list(armaOrder = c(4, 4)),
                      variance.model = list(garchOrder = c(1, 1), model = "sGARCH"),
                      distribution.model = 'norm')
sg11_42 <- ugarchspec(mean.model = list(armaOrder = c(4, 2)),
                      variance.model = list(garchOrder = c(1, 1), model = "sGARCH"),
                      distribution.model = 'norm')
sg11_24 <- ugarchspec(mean.model = list(armaOrder = c(2, 4)),
                      variance.model = list(garchOrder = c(1, 1), model = "sGARCH"),
                      distribution.model = 'norm')
sg11_20 <- ugarchspec(mean.model = list(armaOrder = c(2, 0)),
                      variance.model = list(garchOrder = c(1, 1), model = "sGARCH"),
                      distribution.model = 'norm')
sg11_02 <- ugarchspec(mean.model = list(armaOrder = c(0, 2)),
                      variance.model = list(garchOrder = c(1, 1), model = "sGARCH"),
                      distribution.model = 'norm')

g11_43 <- ugarchfit(sg11_43, rets)
g11_44 <- ugarchfit(sg11_44, rets)
g11_42 <- ugarchfit(sg11_42, rets)
g11_24 <- ugarchfit(sg11_24, rets) # Failure to converge
g11_20 <- ugarchfit(sg11_20, rets)
g11_02 <- ugarchfit(sg11_02, rets)

rm(sg11_43, sg11_44, sg11_42, sg11_24, sg11_20, sg11_02, g11_24)

# Looking at BIC for each model
g11_43 # -6.6799
g11_44 # -6.6794
g11_42 # -6.6819
g11_20 # -6.6859
g11_02 # -6.6859

#   Since the MA(2) and AR(2) have same criteria and the second coefficient is 
# insignificant, we drop it on both models. 

rm(g11_43, g11_44, g11_42, g11_20, g11_02)

sg11_10 <- ugarchspec(mean.model = list(armaOrder = c(1, 0)),
                      variance.model = list(garchOrder = c(1, 1), model = "sGARCH"),
                      distribution.model = 'norm')
sg11_01 <- ugarchspec(mean.model = list(armaOrder = c(0, 1)),
                      variance.model = list(garchOrder = c(1, 1), model = "sGARCH"),
                      distribution.model = 'norm')

g11_10 <- ugarchfit(sg11_10, rets)
g11_01 <- ugarchfit(sg11_01, rets)

rm(sg11_10, sg11_01)

g11_10 # -6.6889 and slightly lower AIC (-6.6995), plus coefficient is significant 
# with Robust Standard Errors at 5%
g11_01 # -6.6889 

rm(g11_01)

# Looking at the residuals
# Histogram
hist(residuals(g11_10), breaks = 50, main ='Histogram of the Residuals for MA(1)-GARCH(1, 1)',
     xlab = NA, )
box()

# QQ plot
qqnorm(residuals(g11_10), ylim = c(-0.05, 0.05)) 
qqline(residuals(g11_10), lwd = 2)
# Clearly not normally distributed

# Changing specification of the model to Student t distributed residuals
sg11_10 <- ugarchspec(mean.model = list(armaOrder = c(1, 0)),
                      variance.model = list(garchOrder = c(1, 1), model = "sGARCH"),
                      distribution.model = 'std')
g11_10 <- ugarchfit(sg11_10, rets)
g11_10 # -6.7599 => better than Normal

rm(sg11_10)
# Looking at the residuals
# Histogram
hist(residuals(g11_10), breaks = 30, main ='Histogram of the residuals')
box()

# QQ plot
qqplot(rt(1000, df = 4.7), as.numeric(residuals(g11_10)), ylab = 'Sample Quantiles', 
       xlab = 'Theoretical Quantiles', main = 'Student\'s t Q-Q Plot')
qqline(as.numeric(residuals(g11_10)))
# Behaves better than Normal distribution

# Testing eGARCH
esg11_10 <- ugarchspec(mean.model = list(armaOrder = c(1, 0)),
                       variance.model = list(garchOrder = c(1, 1), model = "eGARCH"),
                       distribution.model = 'std')
eg11_10 <- ugarchfit(esg11_10, rets)
eg11_10 # -6.7983 => Even better criteria, all highly significant coefficients, 
# non-correlated residuals

rm(esg11_10)

# Looking at the residuals
# Histogram
hist(residuals(eg11_10), breaks = 80, main ='Histogram of the residuals',
     xlim = c(-0.05, 0.05), prob = TRUE)
box()

# QQ plot
qqplot(rt(1000, df = 4.7), as.numeric(residuals(g11_10)), ylab = 'Sample Quantiles', 
       xlab = 'Theoretical Quantiles', main = 'Student\'s t Q-Q Plot')
qqline(as.numeric(residuals(eg11_10)))




