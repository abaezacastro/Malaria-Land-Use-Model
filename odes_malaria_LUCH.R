model_1<-function(t, x, parms, ...){
  
  
  #################################################################
  #define parameters
  #################################################################
  
  phi<-parms[[1]]                #Total production factor
  alpha<-parms[[2]]              #Elasticity land
  beta<-parms[[3]]               #Elasticity capital
  p<-parms[[4]]                  #Crop price
  CA<-parms[[5]]                 #Cost per unit of area planted
  tau<-parms[[6]]                #Cost per case of malaria
  s<-parms[[7]]                  #Income saving rate
  zeta<-parms[[8]]               #Depreciation rate
  f<-parms[[9]]                  #Mosquito fertility
  K0<-parms[[10]]                #carrying capasity per unit of area U
  K0_d<-parms[[10]]*parms[[21]]  #carrying capasity per unit of area D
  HBI<-parms[[11]]               #human biting rate
  cc<-parms[[12]]                #Probability of mosquito being infected from infectious human 
  b<-parms[[13]]                 #Probability of human being infected from infectious mosquito bite
  longevity_M<-parms[[14]]       #days a mosquito live
  delta_M<-1/longevity_M         #Adult mosquito mortality rate
  br_H<-parms[[15]]              #Human population growth rate
  delta_H<-parms[[16]]           #Human mortality rate
  mu<-1/parms[[17]]              #Natural recovery rate
  tau_rate<-parms[[18]]          #Treatment induced recovery rate (back to suceptible)
  sigma<-parms[[24]]             #Loss of inmunity 
  w=parms[[25]]                  #Devilitating effect of malaria on economy
  ds=parms[[26]]                 #Human density U
  ds_d=parms[[27]]               #Human density D
  v=parms[[29]]                  #proportion infected migrant
  l=0                            #rate of land-use change (The land does not change in this model, therefore l=0)
  ro = 0                         #Rate protection againt malaria due to economic improvment (There is not economically portected class in this model, therefore ro and lambda =0)
  lambda = 0                     #Rate loss of protection due to povertt (There is not economically portected class in this model)
  g=parms[[30]]                  #gonotrophyc cycle [days]
  a=HBI / g                      #human biting rate
  #################################################################
  
  #################################################################
  #define state variables
  #################################################################
  C=x[[1]]  #Capital
  X=x[[2]]  #healthy mosquito
  W=x[[3]]  #Infectious mosquito
  S=x[[4]]  #suceptible human
  I=x[[5]]  #infected human
  R=x[[6]]  #inmune (recovered) human
  P=x[[7]]  #Economically protected human
  A=x[[8]]  #Undevelop area
  #################################################################
  
  
  #######################################################################################################################
  #calculate productivity and income
  #######################################################################################################################
  Y = phi *  (A ^ alpha) * (C ^ beta) * ((S + w * I + R + P) ^ (1-alpha-beta))     #Calculate production
  Ic = p * Y - CA * A - tau * I                                              #Calculate Income
  if(Ic<0){Ic=0}
  #######################################################################################################################
  
  
  #######################################################################################################################
  K = K0 * A    #Calculate mosquito carrying capasity in total area U
  Fm = f / g    #Mosqito fertility rate (2.2 is the gonotrophic cycle in days See baeza et al 2015)
  #######################################################################################################################
  

  ##########################################################################################################
  #calculate the differential equations
  ##########################################################################################################
  dA = - l * A                                                                                         #land
  dC = s * Ic - zeta * C                                                                               #capital

  dX = Fm * (X + W) * ( 1 - (X + W) / K ) - a * cc * (I / (I + S + R)) * X                           #mosquito
  dW =  a * cc * (I / (I + S + R)) * X -  W * delta_M                                                #mosquito
  
  dS = br_H  * (S + I + R) + sigma * R + (tau_rate * tau) * I  - a * b * (W / (I + S + R)) * S - delta_H *S  #humans 
  dI = a * b * (W / (I + S + R)) * S  - delta_H * I - (mu+tau_rate * tau) * I                                #humans
  dR = mu * I - delta_H * R  - sigma * R                                                               #humans
  dP = (br_H - delta_H) * P                                                                            #humans
  ###########################################################################################################
  res=c(dC, dX, dW, dS, dI, dR, dP, dA)    #report results of differentials equations
  list(res)
} #function ends here