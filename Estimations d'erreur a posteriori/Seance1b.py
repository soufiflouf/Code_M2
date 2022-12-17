#--------------------------------------------
#   Euler Explicit
#--------------------------------------------

import numpy as np
import math as mt
import matplotlib.pyplot as plt

#   Variables
Dt = [0.0167, 0.0083, 0.0056, 0.0042]  # = [1, 0.5, 0.33, 0.25] seconde
lbd = 1
u0 = 1

errL2 = np.zeros(4)
for i in range(4) :
    
    #   Discretisation du pas de temps
    N= mt.floor(1/Dt[i])
    if (N*Dt[i]!=1):
        Nf=N+1
    else :
        Nf=N
    t=np.zeros(Nf+1)
    for j in range(1,Nf+1):
        t[j]=t[j-1]+Dt[i]
    if (N!=Nf):
        t[Nf]=1
    
    #   Solution exacte et numérique
    sol_ex=np.exp(-lbd*t)*u0
    sol_app=np.zeros(Nf+1)
    sol_app[0]=u0 #Condition initiale u(0)=1
    for j in range(1,Nf+1):
        sol_app[j]=sol_app[j-1]*(1-Dt[i])
    
    #   Erreur L2
    erri=sol_ex-sol_app
    errL2[i]=0
    for j in range(1,Nf):
        errL2[i] = errL2[i] + pow(erri[j],2)
    errL2[i] = mt.sqrt(Dt[i]*errL2[i])
    
plt.plot(np.log(Dt),np.log(errL2),'r')
plt.title('Erreur L2 en fonction du pas de temps',fontsize=10)
plt.plot( np.log(Dt),np.log(Dt),'b')  # Pour verifier que la methode est bien d'ordre 1