#################################################
#TCC nova versao
################################################

##################################################
#exemplo pratico

#funcao
nginar <- function(x,series=NULL)
{
  if (is.null(series))
    series <- deparse(substitute(x))
  xfreq <- frequency(x)
  n <- length(x)
  order.max = 1
  r0 <- acf(x, plot = FALSE)$acf[1]
  r <- acf(x, plot = FALSE)$acf[2:(order.max+1)]
  R <- diag(order.max)
  for(i in 1:order.max){
    for(j in 1:order.max){
      if(i!=j){
        R[i,j] <- r[abs(i-j)]
      }
    }
  }
  residuals <- function(x,coef,mu)
  {
    x<- x
    p<- 1
    e <- NULL
    for(t in (p+1):length(x) )
    {
      e[t] <- x[t] - ( sum(coef*x[(t-1):(t-p)]) + mu)
    }
    return(e)
  }
  coef <- round(solve(R, r), 4)
  xbar <- mean(x)
  mu.e <- xbar*(1-sum(coef)) #mean of the error
  var.error <- r0 - sum(coef*r) #variance of the erro
  mu.x <- mu.e/(1-sum(coef))
  sum.var <- sum( coef*(1-coef))
  Vp <- var.error + mu.x*sum.var
  fitted <- nginar.sim(length(x), alpha = coef, mu = mu.e, n.start=200)
  resid <- residuals(x,coef,mu.e)
  rms <- sqrt(mean(resid^2,na.rm = TRUE))
  AICc. <- n*log(Vp) + n*((1+order.max/n)/(1-(order.max+2)/n))
  AIC. <- n*log(Vp) + 2*order.max
  BIC. <- n*log(Vp) + (order.max/n)*log(n)
  res <- list(order = order.max, coef = coef, mean.e = mu.e, var = var.error, rms = rms,
              fitted.values = fitted, bic= BIC., aicc=AICc., aic = AIC., n.used = n, order.max = order.max,
              resid = resid, method = "yule-walker", series = series, frequency = xfreq, call = match.call())
  class(res) <- "inar"
  return(res)
}
#' Function poinar.sim
#' Simulate from an Inar model
#' @param n the length of outputs series. A strictly positive integer.
#' @param order.max the integer component p is the INAR order.
#' @param alpha a vector of INAR coefficients.
#' @param lambda the mean of the poisson distribution.
#' @param n.start the length of 'burn-in' period. If na, the default, a reasonable valve is
computed.
#'@return A time-series object of class "ts".
#'@seealso \code{\link{poinar}}
#'@references
#' Du, J.G. and Li,Y.(1991): The integer-valued autorregressive (INAR(p)) model.
#' \emph{Journal of time series analysis} \bold{12}, 129--142.
#'@examples
#'# A Poisson INAR simulation
#'ts.sim <- poinar.sim(n = 100, order.max = 2, alpha = c(0.1,0.4),lambda = 2, n.start=200)
#'ts.plot(ts.sim)
#' @export

poinar.sim <- function(n, order.max, alpha,lambda, n.start=NA){
  length. <- n + n.start
  x <- rep(NA, times = length.)
  error <- rpois(length., lambda)
  for (i in 1:order.max) {
    x[i] <- error[i]
  }
  for (t in (order.max + 1):length.) {
    x[t] <- 0
    for (j in 1:order.max) {
      x[t] <- x[t] + rbinom(1, x[t - j], alpha[j])
    }
    x[t] <- x[t] + error[t]
  }
  ts(x[(n.start+1):length.],frequency = 1,start=1)
}

#Apos a leitura da funcao a aplicacao

ts.sim <- poinar.sim(n = 100, order.max = 2, alpha = c(0.1,0.4),lambda = 2, n.start=200)
ts.plot(ts.sim, xlab="Tempo", ylab="Frequência",
        main="Exemplo prático de uso do pacote")

##################################################

######################################################
#Lendo os dados

rm(dados)
dados=read.table(file="dados.csv", header=T, sep=";") 
str(dados)
attach(dados)
names(dados)
######################################################

######################################################
#Analise descritiva

summary(dados)

#grfico da serie
ts.plot(dados$X10_14, xlab="Tempo",ylab="Frequência",main="10|-|14")
#
ts.plot(dados$X15_19, xlab="Tempo",ylab="Frequência",main="15|-|19")
#
ts.plot(dados$X20_24, xlab="Tempo",ylab="Frequência",main="20|-|24")
#
ts.plot(dados$X25_29, xlab="Tempo",ylab="Frequência",main="25|-|29")
#
ts.plot(dados$X30_34, xlab="Tempo",ylab="Frequência",main="30|-|34")
#
ts.plot(dados$X35_39, xlab="Tempo",ylab="Frequência",main="35|-|39")
#
ts.plot(dados$X40_44, xlab="Tempo",ylab="Frequência",main="40|-|44")
#
ts.plot(dados$X45_49, xlab="Tempo",ylab="Frequência",main="45|-|49")
#
ts.plot(dados$Total, xlab="Tempo",ylab="Frequência",main="Total")
#

######################################################
#Gráfico da FAC's e FACP's

par(mfrow=c(2,1))
acf(dados$X10_14,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(dados$X10_14,main="",xlab="Lag",ylab="FACP")
#

##
par(mfrow=c(2,1))
acf(dados$X15_19,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(dados$X15_19,main="",xlab="Lag",ylab="FACP")
#

##
par(mfrow=c(2,1))
acf(dados$X20_24,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(dados$X20_24,main="",xlab="Lag",ylab="FACP")
#

##
par(mfrow=c(2,1))
acf(dados$X25_29,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(dados$X25_29,main="",xlab="Lag",ylab="FACP")
#

##
par(mfrow=c(2,1))
acf(dados$X30_34,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(dados$X30_34,main="",xlab="Lag",ylab="FACP")
#

##
par(mfrow=c(2,1))
acf(dados$X35_39,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(dados$X35_39,main="",xlab="Lag",ylab="FACP")
#

##
par(mfrow=c(2,1))
acf(dados$X40_44,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(dados$X40_44,main="",xlab="Lag",ylab="FACP")
#

##
par(mfrow=c(2,1))
acf(dados$X45_49,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(dados$X45_49,main="",xlab="Lag",ylab="FACP")
#

##
par(mfrow=c(2,1))
acf(dados$Total,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(dados$Total,main="",xlab="Lag",ylab="FACP")
#

#########################################

  