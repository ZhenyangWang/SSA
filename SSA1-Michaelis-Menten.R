## Simple implementation of the Stochastic Simulation Algorithm
## (or Gillespie's algorithm) on the Michaelis-Menten system
## 
## A subsrate S is converted to a product P only in the presence
## of a catalyst E (enzyme).The model for this is:
## S + E -> SE   with stochastic rate k1
## SE -> S + E   with stochastic rate k2
## SE -> S + P   with stochastic rate k3
## 
## This function can plot different simulation by change the copy number of species
## Last updated on 2017-5-19
## Author: ZhenyangWang


MM = function()
{
  ########## Reaction Stoichiometric Matrix ##########
  V0 = c(-1,1,0,-1,1,1,1,-1,-1,0,0,1)                # reaction rules
  V = matrix(V0,byrow=T, nrow=4)
  
  ########## Parameters and Initial Conditions ##########
  X = matrix(0,4,1)
  X[1,] = 200                                        # copy numbers of substrate
  X[2,] = 100                                        # copy numbers of enzyme
  rate = c(1.66e-3,1e-4,0.1)                         # stochastic rates
  m = length(rate)
  a = vector("numeric",m)
  t = 0
  tfinal = 60                                        # reaction times
  
  ########## Initial Reaction State ##########
  S = c()
  E = c()
  SE = c()
  P = c()
  times = c()
  
  S = rbind(S,X[1,])                                 # record the initial number of subsrate
  E = rbind(E,X[2,])                                 # record the initial number of enzyme
  SE = rbind(SE,X[3,])                               # record the initial number of complex
  P = rbind(P,X[4,])                                 # record the initial number of product
  times = rbind(times,0)
  
  ########## While Loop for Reaction Processing ##########
  while (t < tfinal)                                 
  {
    a[1] = rate[1]*X[1,]*X[2,]                       # calculate the propensity for first reaction rule
    a[2] = rate[2]*X[3,]                             # the propensity for the second one
    a[3] = rate[3]*X[3,]
    a0 = sum(a)                                      # the total propensity
    u = sample(1:m,1,replace=FALSE,prob = c(a/a0))   # pick a reaction by calculate probability for firing each reaction, which based on the propensity of each reaction
    X = X+V[,u]                                      # reaction matrix processing
    
    tp = rexp(1,a0)                                  # the time step,which calculate the time of next reaction
    t = t+tp
    
    S = rbind(S,X[1,])                               # record the change of number of subsrate
    E = rbind(E,X[2,])                               # record the change of number of enzyme
    SE = rbind(SE,X[3,])                             # record the change of number of complex
    P = rbind(P,X[4,])                               # record the change of number of product
    times = rbind(times,t)                           # record the time step
  }
  
  ########## Plots ##########
  plot(times,S,ylim=c(min(S,E,SE,P),max(S,E,SE,P)),col='red',type = "l",main = "Michaelis-Menten kinetics",sub = "1X copy number",xlab = "time",ylab = "number of species")
  lines(times,E,col='green')
  lines(times,SE,col='blue')
  lines(times,P,col='cyan')
  legend("topleft", legend=c("S","E","SE","P"),lty=1,col=c("red","green","blue","cyan"),bty = "n")
}
MM()
