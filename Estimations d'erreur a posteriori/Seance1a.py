#-----------------------------------------------------------------------------
#   Euler Explicit
#-----------------------------------------------------------------------------

import numpy as np
import math as mt
import matplotlib.pyplot as plt

#   Variables
Dt = 0.016 # = 1 seconde
lbd = 1
u0 = 1

#   Discretisation du pas de temps
N = mt.floor(1/Dt)
if (N*Dt!=1):
    Nf=N+1
else :
    Nf=N
t=np.zeros(Nf+1)
for i in range(1,Nf+1):
    t[i]=t[i-1]+Dt
if (N!=Nf):
    t[Nf]=1

#   Solution exacte et numérique
sol_ex=np.exp(-lbd*t)*u0
sol_app=np.zeros(Nf+1)
sol_app[0]=u0 #Condition initiale u(0)=1
for j in range(1,Nf+1):
    sol_app[j]=sol_app[j-1]*(1-Dt)

print(sol_app)
print (sol_ex)
    
# #   Plot
plt.plot(t,sol_app,'r',t,sol_ex,'b')
plt.title('Solution exacte et approchée')
