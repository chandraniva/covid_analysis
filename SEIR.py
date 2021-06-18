# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 18:35:05 2021
@author: Chandraniva
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import csv

def model(y,t,beta,N,sigma,gamma,mu):
    [S, I, R, Q] = y
    
    dIdt = beta*I*S/N - gamma*I
    dRdt = gamma*I          
    dSdt = - beta*I*S/N - S*dIdt
    dQdt =  S*dIdt
    dydt = [dSdt, dIdt, dRdt, dQdt]
    
    return dydt

def rk4(x0, y0, x, h,beta,N,sigma,gamma,mu):
    # Count number of iterations using step size or
    # step height h
    n = (int)((x-x0)/h)
    # Iterate for number of iterations
    y = y0
    y_hist=[y0]
    for i in range(1, n+1):
        "Apply Runge Kutta Formulas to find next value of y"
        k1 =  [ h * model(y,x0,beta,N,sigma,gamma,mu)[i] for i in range (4)]
        k2 =  [ h * model( y + 0.5 * k1, x0 + 0.5 * h,beta,N,sigma,gamma,mu)[i] for i in range (4)]
        k3 =  [ h * model( y + 0.5 * k2,x0 + 0.5 * h,beta,N,sigma,gamma,mu)[i] for i in range (4)]
        k4 =  [ h * model( y + k3,x0 + h,beta,N,sigma,gamma,mu)[i] for i in range (4)]
 
        # Update next value of y
        y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        y_hist.append(y)
        
        # Update next value of x
        x0 = x0 + h
    return y_hist

"""
fears = np.linspace(0,1,100)

def fear_fact(I):
    I_cutoff = 1e5
    if I>I_cutoff:
        return 0
    else:
        return 
    
    if I>I_cutoff:
        return ((I-I_cutoff)/(4*1e1*I_cutoff))**(2)+1
    else:
        return 1
    
    
fears = []
In = []

for i in range(0,int(1e7),1000):
    In.append(i)
    fears.append(fear_fact(i))
      
plt.plot(In,fears)
plt.show()
"""

beta= 0.34 #0.34 
sigma= 0.1 
gamma= 0.06
mu= 0.005
N = 1.3526e9
birth_rate, death_rate = 28230/N ,28230/N # 28230/N , 77756/N

y_init = [N, 3, 3, 0]

"""
t = np.linspace(0,1000,500)
sol = odeint(model, y_init, t, full_output=0,args=(beta,N,sigma,gamma,mu,
                                     birth_rate,death_rate))
"""
t = np.linspace(0,500,501)
sol=rk4(0,y_init,500,0.1,beta,N,sigma,gamma,mu)

tot=sol[:,1]+sol[:,0]+sol[:,2]+sol[:,3]
plt.plot(t, sol[:,1], 'b', label='Infectives')
#plt.plot(t,tot, label='total')
plt.plot(t, sol[:,2], 'g', label='Recovered')
plt.plot(t, sol[:,0], 'r', label='Susceptible')
plt.plot(t, sol[:,3], 'black', label='Quarantined')
plt.legend(loc='best')
plt.xlabel('time (days)')
plt.grid()
plt.show()


infectives_obs = []
total_recovered_obs = []
time = []

with open('case_time_series.csv','r') as csvfile:
    t = 0
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        infectives_obs.append(int(row[2]))
        total_recovered_obs.append(int(row[5]))
        time.append(t)
        t += 1
"""        
plt.plot(time[32:],infectives_obs[32:],label="Infectives observed")
#plt.plot(time[32:],total_recovered_obs[32:], label = "Recovered observed")
plt.legend(loc='best')
plt.xlabel('time (days)')
plt.grid()
plt.show()
"""











