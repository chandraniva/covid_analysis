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
        
        dRdt = gamma*I - alpha*(R_hist[-179]-R_hist[-180])  
        dSdt = - beta*I*S/N - S*fear_function(I,dIdt)/N + alpha*(R_hist[-179]-R_hist[-180])   
        
    else:
        dRdt = gamma*I 
        dSdt = - beta*I*S/N - S*fear_function(I,dIdt)/N 
        
    dQdt =  S*fear_function(I, dIdt)/N    
    dydt = [dSdt, dIdt, dRdt, dQdt]
    
    return np.array(dydt)


def rk4(t,y0,t0,I_cutoff):
    h=1e-2
    n = (int)((t - t0)/h)
    y = y0
    R_hist = []
    y_hist = []
    for i in range(1, n+1):
        
        [S, I, R, Q] = y
        
        if i%(int(1/h))==1:
            R_hist.append(R)
            y_hist.append(y)
        
        
        if I<I_cutoff:
            
            S = S+Q
            Q = 0
            y = np.array([S, I, R, Q])
            
            k1 = h * dydt(t0, y, R_hist)
            k2 = h * dydt(t0 + 0.5 * h, y + 0.5 * k1, R_hist)
            k3 = h * dydt(t0 + 0.5 * h, y + 0.5 * k2, R_hist)
            k4 = h * dydt(t0 + h, y + k3, R_hist)
            y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
            
        else:
            
            k1 = h * dydt2(t0, y, R_hist)
            k2 = h * dydt2(t0 + 0.5 * h, y + 0.5 * k1, R_hist)
            k3 = h * dydt2(t0 + 0.5 * h, y + 0.5 * k2, R_hist)
            k4 = h * dydt2(t0 + h, y + k3, R_hist)
            y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        
        t0 = t0 + h
        y = [y[i] for i in range(4)]
        
    return y_hist

alpha = 0.5
beta= 0.08
gamma= 0.02
N = 1.3526e9
I_cutoff = 1e4


S_init, I_init, R_init, Q_init = N-2006, 2000, 6, 0
y_init = [S_init, I_init, R_init, Q_init]

days = 400

sol = np.array(rk4(days,y_init,0,I_cutoff))
print(sol)

plt.figure()
tot=sol[:,1]+sol[:,0]+sol[:,2]+sol[:,3]
plt.plot(np.arange(days), sol[:,1], 'b', label='Infectives')
#plt.plot(np.arange(days),tot, label='total')
#plt.plot(np.arange(days), sol[:,2], 'g', label='Recovered')
#plt.plot(np.arange(days), sol[:,0], 'r', label='Susceptible')
#plt.plot(np.arange(days), sol[:,3], 'black', label='Quarantined')
plt.legend(loc='best')
plt.xlabel('time (days)')
plt.grid()
plt.show()

