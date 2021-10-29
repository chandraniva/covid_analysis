import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import csv

with open('.\data\case_time_series.csv', newline='') as csvfile:
    data = list(csv.reader(csvfile))

active = []
daily_deceased = []

for i in range (1,len(data)):
    active.append(int(data[i][3])-int(data[i][5])-int(data[i][7]))
    daily_deceased.append(int(data[i][6]))

active = np.array(active)
daily_deceased = np.array(daily_deceased)

print(max(active[1:100]))

plt.plot(active[425:464]/max(active[425:518]), daily_deceased[425:464]/max(daily_deceased[425:518]),'x',label="increase")
plt.plot(active[465:518]/max(active[425:518]), daily_deceased[465:518]/max(daily_deceased[425:518]),'x',label="decrease")
plt.xlabel("Active cases")
plt.ylabel("Daily deceased")
plt.legend()
plt.show()

with open('.\data\India_COVID_Mobility.csv', newline='') as csvfile:
    mobility_data = list(csv.reader(csvfile))
    
#print(mobility_data[0])

sum_except_groceries = []

for i in range (1,len(mobility_data)):
    sum_except_groceries.append(float(mobility_data[i][7]))
    
print(data[425][1])    
print(data[518][1])
print(data[464][1])
    
print(mobility_data[409][0]) 
print(mobility_data[502][0]) 
print(mobility_data[448][0]) 

def pol(x,a,p):
    return a*x**p

a_fit, p_fit = curve_fit(pol, daily_deceased[425:518], sum_except_groceries[409:502])[0]
print(a_fit,p_fit)

plt.plot(daily_deceased[425:464],sum_except_groceries[409:448],'x',label="going up")
plt.plot(daily_deceased[465:518],sum_except_groceries[449:502],'x',label="going down")
plt.plot(np.arange(0,6000),pol(np.arange(0,6000),a_fit,p_fit),label="fit")
plt.xlabel("daily dec cases")
plt.ylabel("sum except groceries")
plt.legend()
plt.show()  
 
#print(sum_except_groceries[1][0])
#print(daily_deceased[1][1])

subtracted_sum_ex_gro=[]
    
for i in range(len(sum_except_groceries)):
    subtracted_sum_ex_gro.append(sum_except_groceries[i] - a_fit*daily_deceased[i+16]**p_fit)
    
      
def exp_decay(x,c,d):
    return c*np.exp(-d*x)
    
c_fit, d_fit = curve_fit(exp_decay, np.arange(len(subtracted_sum_ex_gro[90:])), subtracted_sum_ex_gro[90:],p0=[60,1])[0]
print(c_fit,d_fit)    
    
plt.plot(subtracted_sum_ex_gro[90:],'x',label="Subtracted ") 
plt.plot(np.arange(0,410), exp_decay(np.arange(0,410),c_fit,d_fit))
plt.show()    

print(len(subtracted_sum_ex_gro[90:]))

final_subtracted = []
for i in range(len(subtracted_sum_ex_gro[90:])):
    final_subtracted.append(subtracted_sum_ex_gro[90+i]-c_fit*np.exp(-i*d_fit))
    

print(len(final_subtracted))

def sine(x,A,omega,delta):
    return A*np.sin(omega*x+delta)

np.savetxt("wiggles.txt",final_subtracted)

A_fit, omega_fit, delta_fit =  curve_fit(sine, np.arange(len(final_subtracted)), final_subtracted,p0=[6,90,np.pi/2.])[0]
print(A_fit, omega_fit,delta_fit)

plt.plot(final_subtracted,'x')
plt.plot(np.arange(len(final_subtracted)),sine(np.arange(len(final_subtracted)),A_fit,omega_fit,delta_fit))
plt.show()   


    
    
    
    
    
