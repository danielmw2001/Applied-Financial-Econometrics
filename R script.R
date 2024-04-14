#B184777 R-script

#Prelims
library(xts)
library(rugarch)
library(tseries)
library(moments)
library(forecast)
library(gt)
library(fGarch)
library(pracma)

exchange <- read.csv('GBPUSD.csv')
exchange$Date <- as.Date(exchange[,1])
exchange <- exchange[,1:6]
ets <- xts(exchange[, -1], order.by=exchange$Date)
plot(ets$Close, type='l', main='Figure 1: GBP/USD Exchange Rates', xlab='Date', ylab='Close Price')

#B.1
#Generate log returns and returns
ets$lagclose <- lag(ets$Close, 1)
ets$returns <- (ets$Close / ets$lagclose) -1
ets$logreturns <- log(ets$Close / ets$lagclose)

#Skewness/Kurtosis/Departure from normality in returns
skc <- skewness(na.omit(ets$Close))
kuc <- kurtosis(na.omit(ets$Close))
mec <- mean(na.omit(ets$Close))
varc <- var(na.omit(ets$Close))

sk <- skewness(na.omit(ets$returns))
ku <- kurtosis(na.omit(ets$returns))
me <- mean(na.omit(ets$returns))
var <- var(na.omit(ets$returns))

# making table 
headers <- c("Mean", "Variance", "Skewness", "Kurtosis")

data <- data.frame(me, var, sk, ku)
closedata <- data.frame(mec, varc, skc, kuc)
colnames(data) <- headers
colnames(closedata) <- headers

gt_table2 <- gt(data) %>%
  tab_header(
    title = "Fig 3: Return Moment Data",
  )
gt_table <- gt(closedata) %>%
  tab_header(
    title = "Fig 2: Exchange Rate Moment Data",
  )

gt_table
gt_table2

#negative skewness, and leptokurtotic, as expected. Now Jarque-Bera
JB <- length(ets$returns)*((skewness(na.omit(ets$returns))^2/6) +((kurtosis(na.omit(ets$returns))-3)^2)/24)
jb <- jarque.bera.test(na.omit(ets$returns))
#Manual and library find we strongly reject H0 of normality

#Dicky fuller and differencing
adf.test(ets$Close, k = 1)
adf.test(log(ets$Close), k = 1)
#cannot reject H0 of non-stationarity, p val of 0.47
ets$closediff <- diff(ets$Close)
adf.test(na.omit(ets$closediff), k = 1)
#Can reject non stationarity of differenced time series

#ARCH/GARCH
arch <- garchFit(formula = ~ garch(1,0),data = na.omit(ets$logreturns), trace = F)
garch <- garchFit(formula = ~ garch(1,1),data = na.omit(ets$logreturns), trace = F)
#Table
headers <- c("Model", "AIC", "BIC")

data <- data.frame(Model = c("ARCH(1)", "GARCH(1,1)"),
                   AIC = c("-7.49", "-7.58"),
                   BIC = c("-7.47", "-7.57"))

colnames(data) <- headers

gt_table <- gt(data) %>%
  tab_header(
    title = "Figure 4: Return Moment Data"
  )

gt_table

#Hurst Exp long range dependence
hurstexp(na.omit(ets$logreturns)^2)

#B.2
# ADF and differencing on Close series
adf.test(ets$Close, k = 1)
ets$closediff <- diff(ets$Close)
adf.test((na.omit(ets$closediff)), k = 1)

#ACF/PACF
# ACF with title
acf(na.omit(ets$closediff), main = "Figure 5: Autocorrelation Function")

# PACF with title
pacf(na.omit(ets$closediff), main = "Figure 6: Partial Autocorrelation Function")

#Model selection and table generation
current_low_AIC <- Inf  
current_low_BIC <- Inf
current_low_AICc <- Inf
AICs <- list()
AICcs <- list()
BICs <- list()
model <- list()
for (i in 0:8) {
  for (j in 0:8) {
    # Fit the ARIMA model
    arima_fit <- Arima(ets$Close, order = c(i, 1, j))
    if (arima_fit$aic < current_low_AIC) {
      current_low_AIC <- arima_fit$aic
      best_AIC_params <- c(i, 1, j)
    }
    # Check BIC
    if (arima_fit$bic < current_low_BIC) {
      current_low_BIC <- arima_fit$bic
      best_BIC_params <- c(i, 1, j)
    }
    if (arima_fit$aicc < current_low_AICc) {
      current_low_AICc <- arima_fit$aicc
      best_AICc_params <- c(i, 1, j)}
    
    model <- append(model, list(c(i,1,j)))
    BICs <- append(BICs, arima_fit$bic)
    AICs <- append(AICs, arima_fit$aic)
    AICcs <- append(AICcs, arima_fit$aicc)
    
  }
}

library(forecast)

AICs <- c()
BICs <- c()
AICcs <- c()
models <- c()

for (i in 0:2) {
  for (j in 0:2) {
    arima_fit <- Arima(ets$Close, order = c(i, 1, j))
    AICs <- c(AICs, arima_fit$aic)
    BICs <- c(BICs, arima_fit$bic)
    AICcs <- c(AICcs, arima_fit$aicc)
    models <- c(models, paste("ARIMA(", i, ",1,", j, ")", sep = ""))
  }
}

criteria_df <- data.frame(
  Model = models,
  AIC = AICs,
  BIC = BICs,
  AICc = AICcs
)

gt(criteria_df)

#RMSE calculations
training <- ets$Close[1:1013]
test <- ets$Close[1014:1352]

aicarima <- Arima(training, order = c(8,1,2))
bicarima <- Arima(training, order = c(0,1,0))

aicforecast <- forecast(aicarima, h = length(test))
bicforecast <- forecast(bicarima, h = length(test))

aic_forecast_series <- c(training, aicforecast$mean)
bic_forecast_series <- c(training, bicforecast$mean)

bicrmse <- sqrt(mean((bicforecast$mean - as.numeric(test))^2))
aicrmse <- sqrt(mean((aicforecast$mean - as.numeric(test))^2))


true_values <- c(as.numeric(training), as.numeric(test))
time_seq <- seq_along(true_values)

plot(index(ets)[1000:1352], true_values[1000:1352], type = "l", col = "black", xlab = "Time", ylab = "Close", main = "Figure 7: Actual vs. Forecasted Values")

lines(index(ets)[1014:1352], aicforecast$mean, col = "blue", lty = 2)

lines(index(ets)[1014:1352], bicforecast$mean, col = "red", lty = 2)

legend("topleft", legend = c("Actual", "AIC Forecast", "BIC Forecast"), col = c("black", "blue", "red"), lty = 1:2, cex = 0.8)

#Compare if AIC model RMSE is lower for short windows
counter <- 0
counter2 <- 0
for (i in 0:5){
  training <- ets$Close[(1+(5*i)):(1013+(5*i))]
  test <- ets$Close[(1014+(5*i)):(1019+(5*i))]
  
  aicarima <- Arima(training, order = c(8,1,2))
  bicarima <- Arima(training, order = c(0,1,0))
  aicforecast <- forecast(aicarima, h = length(test))
  bicforecast <- forecast(bicarima, h = length(test))
  for (j in 1:5){
    aicshortrmse <- sqrt(mean((aicforecast$mean[1:j] - as.numeric(test[1:j]))^2))
    bicshortrmse <- sqrt(mean((bicforecast$mean[1:j] - as.numeric(test[1:j]))^2)) 
    counter2 <- counter2 +1
    if (aicshortrmse < bicshortrmse){
      counter <- counter + 1}
  }
}

#Final Forecasting
final <- Arima(ets$Close[1000:1352], order = c(0,1,0))
forecast <- forecast(final, h = 5)  
plot(forecast, main = "Figure 6: Forecast GBP/USD Exchange Rates", ylab = "Close Price")

#B.3
#Historical sim Variance method
VaR_99 <- quantile(na.omit(ets$logreturns), 0.01)
sample_mean <- mean(na.omit(ets$returns))
sample_sd <- sd(na.omit(ets$returns))
z_score <- qnorm(0.01)
VaR_1_percent <- sample_mean + z_score * sample_sd

#Monte Carlo 1 day and 5 day simulation
returns <- ets$returns
fx <- auto.arima(ets$returns, approximation = F)
spread <- list()
fivedayspread <- list()
fc <- forecast(fx, h = 1)

mean_forecast <- as.numeric(fc$mean)

for (i in (1:1000)){
  spread <- append(spread, (mean_forecast + sample(fx$residuals, size = 1)))}

for (i in 1:1000) {
  for (j in 1:5) {
    x <- Arima(returns, order = c(0,0,0))
    fc_value <- forecast(x, h = 1)$mean
    sampled_residual <- sample(fx$residuals, size = 1)
    last_date <- last(index(returns))
    next_date <- seq(last_date, by = "days", length.out = 2)[2]
    
    new_xts <- xts(fc_value + sampled_residual, order.by = next_date)
    
    
    returns <- rbind(returns, new_xts)
  }
  sum_returns <- sum(as.numeric(returns[(length(returns)-4):length(returns)]))
  fivedayspread <- append(fivedayspread, sum_returns)
  returns <- head(returns, -5)
  print(i)
}

VaR <- quantile(unlist(spread), 0.01)
fivedayVaR <- quantile(unlist(fivedayspread), 0.01)

#B.4
#Dummy creation and model fitting
ets$prepandemic <- as.integer(format(index(ets), "%Y-%m-%d") >= "2020-03-04" & format(index(ets), "%Y-%m-%d") < "2020-03-11")
ets$pandemicday <- as.integer(format(index(ets), "%Y-%m-%d") =="2020-03-11")
ets$pandemic <- as.integer(format(index(ets), "%Y-%m-%d") >="2020-03-12" & format(index(ets), "%Y-%m-%d") < "2020-03-18")
ets$postpandemic <- as.integer(format(index(ets), "%Y-%m-%d") >="2020-03-18" & format(index(ets), "%Y-%m-%d") < "2020-03-25")
ets$pandemicover <- as.integer(format(index(ets), "%Y-%m-%d") =="2023-05-05")
covid <- auto.arima(ets$returns, approximation = F, xreg = cbind(ets$pandemic, ets$prepandemic, ets$postpandemic, ets$pandemicday,ets$pandemicover))

#Abnormal Returns Calculations
expectedmodel <- auto.arima(na.omit(ets$returns["/2020-03-10"]), approximation = F)
expectedreturns <- forecast(expectedmodel, 10)$mean
actualreturns <- ets$logreturns["2020-03-11/2020-03-24"]
abnormalreturns <- as.numeric(actualreturns) -as.numeric(expectedreturns) 
mean(abnormalreturns)
t_stat <- t.test(abnormalreturns)
t_stat
