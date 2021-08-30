import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def fear_function(I,dIdt,t):
    beta2 = 0.1 
    fear = I**2*beta2 - t*100
    return fear


def dydt(t,y):
    [S, I, R, Q] = y
    dIdt = beta*I*S/N - gamma*I
    dRdt = gamma*I          
    dSdt = - beta*I*S/N 
    dQdt = 0    
    dydt = [dSdt, dIdt, dRdt, dQdt]
    
    return np.array(dydt)

def dydt2(t,y):
    [S, I, R, Q] = y
    dIdt = beta*I*S/N - gamma*I
    dRdt = gamma*I          
    dSdt = - beta*I*S/N - S*fear_function(I,dIdt,t)/N
    dQdt =  S*fear_function(I, dIdt,t)/N    
    dydt = [dSdt, dIdt, dRdt, dQdt]
    
    return np.array(dydt)


def rk4(t,y0,t0,I_cutoff):
    h=1e-1
    n = (int)((t - t0)/h)
    y = y0
    for i in range(1, n+1):
        [S, I, R, Q] = y
        
        if I<I_cutoff:
            #print("#",end="")
            S = S+Q
            Q = 0
            y = np.array([S, I, R, Q])
            
            k1 = h * dydt(t0, y)
            k2 = h * dydt(t0 + 0.5 * h, y + 0.5 * k1)
            k3 = h * dydt(t0 + 0.5 * h, y + 0.5 * k2)
            k4 = h * dydt(t0 + h, y + k3)
            y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
            
        else:
            #print(".",end="")
            k1 = h * dydt2(t0, y)
            k2 = h * dydt2(t0 + 0.5 * h, y + 0.5 * k1)
            k3 = h * dydt2(t0 + 0.5 * h, y + 0.5 * k2)
            k4 = h * dydt2(t0 + h, y + k3)
            y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        
        t0 = t0 + h
        
    return y

beta= 0.34  
gamma= 0.1
N = 1.3526e9
I_cutoff = 1e4

y_init = [N-2006, 2000, 6, 0]

days = 300
sol = np.zeros((days,4))
t = np.zeros(days)


for i in range (days):
    sol[i,:] = rk4(i,y_init,0,I_cutoff)
    t[i] = i

print(sol)

plt.figure()
tot=sol[:,1]+sol[:,0]+sol[:,2]+sol[:,3]
plt.plot(t, sol[:,1], 'b', label='Infectives') 
#plt.plot(t,tot, label='total')
#plt.plot(t, sol[:,2], 'g', label='Recovered')
#plt.plot(t, sol[:,0], 'r', label='Susceptible')
#plt.plot(t, sol[:,3], 'black', label='Quarantined')
plt.legend(loc='best')
plt.xlabel('time (days)')
plt.grid()
plt.show()