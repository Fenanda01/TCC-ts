#################################################
#TCC nova versao
################################################
#pacotes para utilizar o GIT e/ou Git\hub
install.packages("usethis")
library(usethis)
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

#######################
dadossup=dados
#######################
names(dados)

Ano=dados$Ano
Mes=dados$Mes
##########################

#grafico da serie

###############################################################
ts.plot(dados$X10_14, xlab="Tempo",ylab="Frequência",main="10|-|14")
#Nao estacionario
acf(dados$X10_14,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#O decaimento e linear e lento, dando seguranca para dizer que a serie e nao estacionaria

#Armazendo os dados em uma variavél de suporte para aplicar a diferenca
#Aplicando o filtro primeira vez (diferenca)
A=diff(dados$X10_14)
dados["Classe1"]=A
ts.plot(A, xlab="Tempo",ylab="Frequência",main="Gráfico da 1ª diferença da série 10 a 14")
#Com a primeira diferenca se tornou estacionario

###############################################################
ts.plot(dados$X15_19, xlab="Tempo",ylab="Frequência",main="15|-|19")
#Nao estacionario
acf(dados$X15_19,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#O decaimento e linear e lento, dando seguranca para dizer que a serie e nao estacionaria

#Armazendo os dados em uma variavél de suporte para aplicar a diferenca
#Aplicando o filtro primeira vez (diferenca)
B=diff(dados$X15_19)
dados["Classe2"]=B
ts.plot(B, xlab="Tempo",ylab="Frequência",main="Gráfico da 1ª diferença da série 15 a 19")
#Com a primeira diferenca se tornou estacionario

###############################################################
ts.plot(dados$X20_24, xlab="Tempo",ylab="Frequência",main="20|-|24")
#Nao estacionario
acf(dados$X20_24,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#O decaimento e linear e lento, dando seguranca para dizer que a serie e nao estacionaria

#Armazendo os dados em uma variavél de suporte para aplicar a diferenca
#Aplicando o filtro primeira vez (diferenca)
C=diff(dados$X20_24)
dados["Classe3"]=C
ts.plot(C, xlab="Tempo",ylab="Frequência",main="Gráfico da 1ª diferença da série 20 a 24")
#Com a primeira diferenca se tornou estacionario

###############################################################
ts.plot(dados$X25_29, xlab="Tempo",ylab="Frequência",main="25|-|29")
#Nao estacionario
acf(dados$X25_29,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#O decaimento e linear e lento, dando seguranca para dizer que a serie e nao estacionaria

#Armazendo os dados em uma variavél de suporte para aplicar a diferenca
#Aplicando o filtro primeira vez (diferenca)
D=diff(dados$X25_29)
dados["Classe4"]=D
ts.plot(D, xlab="Tempo",ylab="Frequência",main="Gráfico da 1ª diferença da série 25 a 29")
#Com a primeira diferenca se tornou estacionario

###############################################################
ts.plot(dados$X30_34, xlab="Tempo",ylab="Frequência",main="30|-|34")
#Nao estacionario
acf(dados$X30_34,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#O decaimento e linear e lento, dando seguranca para dizer que a serie e nao estacionaria

#Armazendo os dados em uma variavél de suporte para aplicar a diferenca
#Aplicando o filtro primeira vez (diferenca)
E=diff(dados$X30_34)
dados["Classe5"]=E
ts.plot(E, xlab="Tempo",ylab="Frequência",main="Gráfico da 1ª diferença da série 30 a 34")
#Com a primeira diferenca se tornou estacionario

###############################################################
ts.plot(dados$X35_39, xlab="Tempo",ylab="Frequência",main="35|-|39")
#Nao estacionario
acf(dados$X35_39,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#O decaimento e linear e lento, dando seguranca para dizer que a serie e nao estacionaria

#Armazendo os dados em uma variavél de suporte para aplicar a diferenca
#Aplicando o filtro primeira vez (diferenca)
F=diff(dados$X35_39)
dados["Classe6"]=F
ts.plot(F, xlab="Tempo",ylab="Frequência",main="Gráfico da 1ª diferença da série 35 a 39")
#Com a primeira diferenca se tornou estacionario

###############################################################
ts.plot(dados$X40_44, xlab="Tempo",ylab="Frequência",main="40|-|44")
#Nao estacionario
acf(dados$X40_44,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#O decaimento e linear e lento, dando seguranca para dizer que a serie e nao estacionaria

#Armazendo os dados em uma variavél de suporte para aplicar a diferenca
#Aplicando o filtro primeira vez (diferenca)
G=diff(dados$X40_44)
dados["Classe7"]=G
ts.plot(F, xlab="Tempo",ylab="Frequência",main="Gráfico da 1ª diferença da série 40 a 44")
#Com a primeira diferenca se tornou estacionario

###############################################################
ts.plot(dados$X45_49, xlab="Tempo",ylab="Frequência",main="45|-|49")
H=dados$X45_49
#estacionario
acf(dados$X45_49,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))

###############################################################
ts.plot(dados$Total, xlab="Tempo",ylab="Frequência",main="Total")
#Nao estacionario
acf(dados$Total,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#O decaimento e linear e lento, dando seguranca para dizer que a serie e nao estacionaria

#Armazendo os dados em uma variavél de suporte para aplicar a diferenca
#Aplicando o filtro primeira vez (diferenca)
I=diff(dados$Total)
dados["Classe9"]=I
ts.plot(I, xlab="Tempo",ylab="Frequência",main="Gráfico da 1ª diferença da série Total")
#Com a primeira diferenca se tornou estacionario

#######################################################

#Apos aplicar a primeira diferenca todas as series se tornaram estacionarias

######################################################
#Os graficos acf devem ser feitos da diferenca 

#Gráfico da FAC's e FACP's
#A ordem do MA eu vejo na FAC
#A ordem do RA eu vejo na FACP

##############################################################
#Identificar se há componentes auto regressivos no modelo
#Nesse passo apenas diz como analisar o gráfico da FAC
#Lembrando que a analise deve ser feita a partir  do lag 1
#para o modelo ser AR (auto regressivo) deve atender a todos os seuintes preceitos
#O decaimento das FAC's é exponencial;
#O tipo de decaimento depende do sinal de ??(fi);
#Quando ?? < 0 as FAC's alternam entre valores positivos e negativos
#Quando ?? > 0 as FAC's são sempre > 0;
#Quanto menor o valor absoluto de ?? (|??|), mais rapidamente as FAC's irão para zero.
##############################################################

######################################################
#CLASSE 1: 10 A 14
par(mfrow=c(2,1))
acf(A,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(A,main="",xlab="Lag",ylab="FACP")
#

######################################################
#CLASSE 2: 15 A 19
##
par(mfrow=c(2,1))
acf(B,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(B,main="",xlab="Lag",ylab="FACP")
#

######################################################
#CLASSE 3: 20 A 24 
##
par(mfrow=c(2,1))
acf(C,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(C,main="",xlab="Lag",ylab="FACP")
#

######################################################
#CLASSE 4: 25 A 29
##
par(mfrow=c(2,1))
acf(D,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(D,main="",xlab="Lag",ylab="FACP")
#

######################################################
#CLASSE 5: 30 A 34
##
par(mfrow=c(2,1))
acf(E,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(E,main="",xlab="Lag",ylab="FACP")
#

######################################################
#CLASSE 6: 35 A 39
##
par(mfrow=c(2,1))
acf(F,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(F,main="",xlab="Lag",ylab="FACP")
#

######################################################
#CLASSE 7 40 A 44
##
par(mfrow=c(2,1))
acf(G,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(G,main="",xlab="Lag",ylab="FACP")
#

######################################################
#CLASSE 8 45 A 49 
##
par(mfrow=c(2,1))
acf(H,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(H,main="",xlab="Lag",ylab="FACP")
#

######################################################
#CLASSE 9: TOTAL
##
par(mfrow=c(2,1))
acf(I,main="",xlab="Lag",ylab="FAC",xlim=c(1,20))
#
pacf(I,main="",xlab="Lag",ylab="FACP")
#

######################################################
summary(A)
######################################################
summary(B)
######################################################
summary(C)
######################################################
summary(D)
######################################################
summary(E)
######################################################
summary(F)
######################################################
summary(G)
######################################################
summary(H)
######################################################
summary(I)
######################################################


ts.sim <- poinar.sim(n = 100, order.max = 2, alpha = c(0.1,0.4),lambda = 2, n.start=200)
ts.plot(ts.sim, xlab="Tempo", ylab="Frequência",
        main="Exemplo prático de uso do pacote")
######################################################
  