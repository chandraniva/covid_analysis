import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def fear_function(I,dIdt):
    beta2 = 0.001 
    fear = beta2*(I-I_cutoff)**2   # dIdt**2
    return fear


def dydt(t,y,R_hist):
    
    [S, I, R, Q] = y
    
    dSdt, dRdt = 0., 0.
    dIdt = beta*I*S/N - gamma*I
    
    if len(R_hist)>180:
        #print(".",end="")
        dRdt = gamma*I - alpha*(R_hist[-179]-R_hist[-180])    
        dSdt = - beta*I*S/N + alpha*(R_hist[-179]-R_hist[-180])    
    
    else:
        dRdt = gamma*I
        dSdt = - beta*I*S/N 
    
    dQdt = 0    
    dydt = [dSdt, dIdt, dRdt, dQdt]
    
    return np.array(dydt)

def dydt2(t,y,R_hist):
    
    [S, I, R, Q] = y
    
    dSdt, dRdt = 0., 0.
    dIdt = beta*I*S/N - gamma*I
    
    if len(R_hist)>180:
        #print("#",end="")
        dRdt = gamma*I - alpha*(R_hist[-179]-R_hist[-180])  
        dSdt = - beta*I*S/N - S*fear_function(I,dIdt)/N + alpha*(R_hist[-179]-R_hist[-180])   
        
    else:
        dRdt = gamma*I 
        dSdt = - beta*I*S/N - S*fear_function(I,dIdt)/N 
        
    dQdt =  S*fear_function(I, dIdt)/N    
    dydt = [dSdt, dIdt, dRdt, dQdt]
    
    return np.array(dydt)


def rk4(t,y0,t0,I_cutoff):
    h=1e-1
    n = (int)((t - t0)/h)
    y = y0
    R_hist = []
    for i in range(1, n+1):
        
        [S, I, R, Q] = y
        R_hist.append(R)
        
        if I<I_cutoff:
            #print("#",end="")
            S = S+Q
            Q = 0
            y = np.array([S, I, R, Q])
            
            k1 = h * dydt(t0, y, R_hist)
            k2 = h * dydt(t0 + 0.5 * h, y + 0.5 * k1, R_hist)
            k3 = h * dydt(t0 + 0.5 * h, y + 0.5 * k2, R_hist)
            k4 = h * dydt(t0 + h, y + k3, R_hist)
            y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
            
        else:
            #print(".",end="")
            k1 = h * dydt2(t0, y, R_hist)
            k2 = h * dydt2(t0 + 0.5 * h, y + 0.5 * k1, R_hist)
            k3 = h * dydt2(t0 + 0.5 * h, y + 0.5 * k2, R_hist)
            k4 = h * dydt2(t0 + h, y + k3, R_hist)
            y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        
        t0 = t0 + h
        
    return y

alpha = 100
beta= 0.34  
gamma= 0.1
N = 1.3526e9
I_cutoff = 1e4


S_init, I_init, R_init, Q_init = N-2006, 2000, 6, 0
y_init = [S_init, I_init, R_init, Q_init]

days = 100
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

