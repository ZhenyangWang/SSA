## Simple implementation of the Stochastic Simulation Algorithm
## (or Gillespie's algorithm) on the Riboflavin pathway model
## 
## The subsrate R5P is converted to a product Riboflavin in the presence
## of several enzyme(RibA,RibH,RibE).The model for this is:
## R5P + RibA <-> R5P.RibA                    with stochastic rate k1,k2
## R5P.RibA -> D2B4P + RibA                   with stochastic rate k3
## D2B4P + RibH <-> D2B4P.RibH                with stochastic rate k1,k2
## D2B4P.RibH-> D8RL + RibH                   with stochastic rate k3
## D8RL + RibE <-> D8RL.RibE                  with stochastic rate k1,k2
## D8RL.RibE -> Riboflavin + RibE             with stochastic rate k3
## 
## This function can plot different simulation by change the copy number of species
## Last updated on 2017-5-19
## Author: ZhenyangWANG.

Ribo = function()
{
  ########## Reaction Stoichiometric Matrix ##########
  V0 = c(-1,1,0,0,0,0,0,0,0,                         # reaction rules
         -1,1,1,0,0,0,0,0,0,
         1,-1,-1,0,0,0,0,0,0,
         0,0,1,-1,1,0,0,0,0,
         0,0,0,-1,1,1,0,0,0,
         0,0,0,1,-1,-1,0,0,0,
         0,0,0,0,0,1,-1,1,0,
         0,0,0,0,0,0,-1,1,1,
         0,0,0,0,0,0,1,-1,-1,
         0,0,0,0,0,0,0,0,1
         )
  V = matrix(V0,byrow=T, nrow=10)
  
  ########## Parameters and Initial Conditions ##########
  X = matrix(0,10,1)
  X[1,] = 301                                        # copy numbers of R5P
  X[2,] = 120                                        # copy numbers of RibA
  X[5,] = 120                                        # copy numbers of RibH
  X[8,] = 120                                        # copy numbers of RibE
  rate = c(1.66e-3,1e-4,0.1)                         # reaction rates
  m = 3*length(rate)
  a = vector("numeric",m)
  t = 0
  tfinal = 100                                       # reaction times
  
  ########## Initial Reaction State ##########
  R5P = c()
  D2B4P = c()
  D8RL = c()
  Riboflavin = c()
  times = c()
  
  R5P = rbind(R5P,X[1,])                             # Record the initial number of R5P
  D2B4P = rbind(D2B4P,X[4,])                         # Record the initial number of D2B4P
  D8RL = rbind(D8RL,X[7,])                           # Record the initial number of D8RL
  Riboflavin = rbind(Riboflavin,X[10,])              # Record the initial number of Riboflavin
  times = rbind(times,0)
  
  ########## While Loop for Reaction Processing ##########
  while (t < tfinal)                                 # reaction keep processing until the set time
  {
    a[1] = rate[1]*X[1,]*X[2,]                       # calculate the propensity for first reaction rule
    a[2] = rate[2]*X[3,]                             # the propensity for the second one
    a[3] = rate[3]*X[3,]                             # the propensity for the third one
    a[4] = rate[1]*X[4,]*X[5,]                       # *
    a[5] = rate[2]*X[6,]                             # *
    a[6] = rate[3]*X[6,]                             # *
    a[7] = rate[1]*X[7,]*X[8,]                       # *
    a[8] = rate[2]*X[9,]                             # *
    a[9] = rate[3]*X[9,]                             # the propensity for the ninth one
    a0 = sum(a)                                      # the total propensity
    u = sample(1:m,1,replace=FALSE,prob = c(a/a0))   # pick a reaction by calculate probability for firing each reaction, which based on the propensity of each reaction
    X = X+V[,u]                                      # reaction matrix processing
    
    tp = rexp(1,a0)                                  # the time step,which calculate the time of next reaction
    t = t+tp
    
    R5P = rbind(R5P,X[1,])                           # record the change of number of R5P
    D2B4P = rbind(D2B4P,X[4,])                       # record the change of number of D2B4P
    D8RL = rbind(D8RL,X[7,])                         # record the change of number of D8RL
    Riboflavin = rbind(Riboflavin,X[10,])            # record the change of number of Riboflavin
    times = rbind(times,t)                           # record the time step
  }
  
  ########## Plots ##########
  plot(times,R5P,ylim=c(min(R5P,D2B4P,D8RL,Riboflavin),max(R5P,D2B4P,D8RL,Riboflavin)),col='red',type = "l",main = "Riboflavin Pathway Model",sub = "1X copy number",xlab = "time",ylab = "number of species")
  lines(times,D2B4P,col='green')
  lines(times,D8RL,col='blue')
  lines(times,Riboflavin,col='cyan')
  legend("topleft", legend=c("R5P","D2B4P","D8RL","Riboflavin"),lty=1,col=c("red","green","blue","cyan"),bty = "n")
}
Ribo()
