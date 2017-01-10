model_2<-function(t, x, parms, ...){
  
  #read parameters
  phi<-parms[[1]]              #total factor of productivity U
  phi_d<-phi*parms[[22]]       #total factor of productivity D
  alpha<-parms[[2]]            #elasticity land
  beta<-parms[[3]]             #elasticity Capital
  p<-parms[[4]]                #crop price
  CA<-parms[[5]]               #cost per unit of area planted
  tau<-parms[[6]]              #cost per case of malaria
  s<-parms[[7]]                #saving rate
  zeta<-parms[[8]]             #depreciation rate
  f<-parms[[9]]                #mosquito fertility
  K0<-parms[[10]]              #carrying capasity per unit of area
  K0_d<-parms[[10]]*parms[[21]]#carrying capasity per unit of area
  g=parms[[30]]
  HBI<-parms[[11]]               
  a=HBI/g                      #human biting rate
  aa_d<-a*parms[[21]]          #human biting rate
  cc<-parms[[12]]              #probability H->M
  b<-parms[[13]]               #probability M->H
  
  LONGEVITY_m<-parms[[14]]     #Moquito lifetime [days]
  delta_M <-1/LONGEVITY_m      #Adult mosquito mortality rate
  
  br_H<-parms[[15]]            #human population growth rate
  delta_H<-parms[[16]]         #human mortality rate
  mu<-1 / parms[[17]]              #natural recovery rate (to R)
  tau_rate<-parms[[18]]        #treatment induced recovery rate (back to suceptible)
  r0<-parms[[19]]              #rate of people moving from S and I to economically protected class  (R) per unit of capital
  p_0<-parms[[20]]             #rate of people moving from R to S per unit of capital
  m<-parms[[23]]               #human movement between regions
  sigma<-parms[[24]]           #lost of inmunity rate
  w=parms[[25]]                #devilitating effect of disease on the economy
  ds=parms[[26]]               #human density U
  ds_d=parms[[27]]             #human density D
  v=parms[[29]]                #proportion infected migrant
  l=parms[[28]]  #rate of land-use change
  
  #define state variables
  
  C=max(c(0,x[[1]]))       #capital U
  X=x[[2]]       #Adult healthy mosquito U
  W=x[[3]]       #Adult infectious mosquito U
  S=x[[4]]       #Suceptible humans U
  I=x[[5]]       #Infected humans U 
  R=x[[6]]       #Inmune humans U
  P=x[[7]]       #Economically protected population U
  
  C_d=max(c(0,x[[8]]))     #Capital D
  X_d=x[[9]]     #Adult healthy mosquito D
  W_d=x[[10]]    #Adult infectious mosquito D
  S_d=x[[11]]    #Suceptible humans D
  I_d=x[[12]]    #Infected humans D
  R_d=x[[13]]    #Inmune humans D
  P_d=x[[14]]    #Economically protected population D
  A=x[[15]]      #area undeveloped
  A_d=x[[16]]    #area developed
  
  #######################################################################################################################
  
  
  
  #######################################################################################################################
  #calculate productivity and income
  Y = phi *  ((A ^ alpha) * (C ^ beta) * ((S + w * I + R + P) ^ (1-alpha-beta)))                        #production in undeveloped
  Y_d = phi_d * ((A_d ^ alpha) * (C_d ^ beta) * ((S_d + w * I_d + R_d + P_d) ^ (1-alpha-beta)))          #production in undeveloped  
  
  Ic = p * Y - CA * A - tau * I          #Income U
  Ic_d = p * Y_d - CA * A_d - tau * I_d  #Income D
  if(Ic<0){Ic=0}                        #condition in case negative income
  if(Ic_d<0){Ic_d=0}                    #condition in case negative income
  #######################################################################################################################
  
  
  
  #######################################################################################################################
  
  ro = 0
  lambda = 0
  
  N_d<-(S_d+I_d+R_d+P_d)
  if(N_d==0){ #if develop area = 0 then there is not economically protected people
    ro_d=0
    lambda_d=0
  }
  else{
    ro_d = r0 * C_d/(S_d+I_d+R_d+P_d)        #rate people accire protected due to impruvement of socio-economic conditions (P) 
    lambda_d = p_0 * (S_d+I_d+R_d+P_d)/C_d   #rate people lost protection due to poverty (e.g. low capital)
  }
  #######################################################################################################################    
  K = K0 * A                               #carrying capasity of total area U
  K_d = K0_d * A_d                         #carrying capasity of total area D
  Fm = f / g                               #mosqito fertility rate (2.2 is the gonotrophic cycle in days See baeza et al 2015)
  #
  #######################################################################################################################
  
  
  #######################################################################################################################
  #differential equation in undeveloped region
  #######################################################################################################################
  
  dA = - l * A                      #land-use change
  
  dC = s * Ic - zeta * C            # capital   
  
  dX = Fm * (X + W) * ( 1 - (X + W) / K ) - a * cc * (I / (I + S + R + P)) * X - l * A * K0 * (X / (X + W))   # healthy mosquito
  dW =  a * cc * (I / (I + S + R + P)) * X - W * delta_M                                                      # infectious mosquito
  
  dS = br_H  * (S + I + R) + sigma * R + lambda * P + (tau_rate * tau) * I  - a * b * (W / (I + S + R + P)) * S - delta_H * S - ro * S  - m * S + m * S_d - l * A * ds * (S/(S+I+R+P)) #human
  dI = a * b * (W / (I + S + R + P)) * S  - delta_H * I - (mu + tau_rate * tau) * I - m * I + m * I_d - l * A * ds * (I/(S+I+R+P))                                                               #human
  dR = mu * I - delta_H * R + m * R_d - m * R - ro * R - sigma * R - l * A * ds * (R/(S+I+R+P))                                                                                        #human
  dP = (br_H - delta_H) * P + ro * (S + R) - lambda * P  - m * P - l * A * ds * (P/(S+I+R+P))                                                                                          #human
  #######################################################################################################################

  #######################################################################################################################    
  #differential equation in developed region
  #######################################################################################################################      
  
  dA_d = l * A
  
  if(A_d>0){  # this equations only change when land is being tranformed (A_d>0) and when people are living in this area (N_d>0)
    if(N_d>0){
      dC_d = s * Ic_d - zeta * C_d                             #capital in developed region
      fi_M = aa_d * cc * (I_d / (I_d + S_d + R_d + P_d))       # force of infection in developed region

      }
    else{
      dC_d = 0   #capital in developed region 
      fi_M=0     # force of infection in developed region
    }
    
    dX_d = Fm * (X_d + W_d) * ( 1 - (X_d + W_d) / K_d ) -  fi_M * X_d + l * A * K0 * (X / (X + W))  #healthy mosquito
    dW_d =  fi_M * X_d - W_d *delta_M  + l * A * K0 * (W / (X + W))                                 #Infectious mosquito 
  }
  else{
    dX_d =0
    dW_d =0
    dC_d = 0 
  }
  
  
  #Calculate force of infection  
  M_d=X_d+W_d  #total mosquito population in developed region
  if(M_d==0){    
    fi_H = 0    #when no mosquitos no tranmission (e.g. force of infection =0)
  }
  else{
    fi_H=aa_d * b * W_d / (I_d + S_d + R_d + P_d)  #Calculate force of infection  
  }
  
  if(A_d>0){ #this statement controls changes in diferential equation only when land started to be tranformed
    dS_d = br_H  * (S_d + I_d + R_d) + sigma * R_d + lambda_d * P_d  +  (tau_rate * tau) * I_d - fi_H * S_d - delta_H * S_d  - ro_d * S_d   - m * S_d + m * S + l * A * ds * (S/(S+I+R+P)) + l * A * (ds_d-ds) * v   #human  
    dI_d = fi_H  * S_d - delta_H * I_d  -  (mu + tau_rate * tau) * I_d - m * I_d + m * I + l * A * ds * (I/(S+I+R+P))                                                                                                #human
    dR_d =  mu * I_d -  delta_H * R_d - m * R_d + m * R - ro_d * R_d + l * A * ds * (R/(S+I+R+P))  - sigma * R_d + (1 - v) * l * A * (ds_d-ds)                                                                       #human
    dP_d =  (br_H - delta_H) * P_d + ro_d * (S_d + R_d) - lambda_d * P_d  + l * A * ds * (P/(S+I+R+P))                                                                                              #human
  }
  else{ #otherwise not changes over time 
    dS_d = 0
    dI_d = 0
    dR_d = 0
    dP_d = 0
  }
  #######################################################################################################################
  
  #report a list containing the differentials 
  res=c(dC, dX, dW, dS, dI, dR, dP, dC_d, dX_d, dW_d, dS_d, dI_d, dR_d, dP_d, dA, dA_d) 
  list(res)
}  #function ends here