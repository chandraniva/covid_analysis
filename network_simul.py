import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si
import matplotlib.animation as animation
import matplotlib.cm as cm
print("\n-------------------------- RHO --------------------------------------\n")

num = 11
N = 1000
tot_population = N*np.ones((num,num))

def gauss_2d(x,y,x0,y0,sigma,c):
    return np.exp(-1/(2*sigma**2)*((x-x0)**2+(y-y0)**2))+c

X, Y = np.meshgrid(np.arange(num), np.arange(num))
rho = gauss_2d(X,Y,(num-1)/2.,(num-1)/2.,3,0.1)
plt.imshow(rho)
plt.colorbar()
plt.title("rho")
plt.show()
  
print("\n-------------------------- evolving infectives ------------------------------------\n")

def sum_neighbours(A):
    n = len(A)
    S = np.zeros_like(A)
    for i in range (n):
        
        k_i, l_i = i-1, i+1
        if i == 0:
            k_i = n-1
        elif i == n-1:
            l_i = 0
                    
        for j in range (n):
            
            k_j, l_j = j-1, j+1
            if j == 0:
                k_j = n-1
            elif j == n-1:
                l_j = 0

            S[i,j] = A[l_i,j]+A[k_i,j]+A[i,l_j]+A[i,k_j]+A[l_i,l_j]+A[l_i,k_j]+A[k_i,l_j]+A[k_i,k_j]
            
    return S
  

def network(y,t,beta,gamma):
    S = y[0:num**2].reshape((num,num))
    I = y[num**2:2*num**2].reshape((num,num))
    R = y[2*num**2:3*num**2].reshape((num,num))
    
    dIdt = beta*rho*S*I/N - gamma*I + beta*rho*S*sum_neighbours(I)/N 
    dRdt = gamma * I
    dSdt = -beta*rho*S*I/N - beta*rho*S*sum_neighbours(I)/N 
    
    dydt = np.concatenate((dSdt.reshape(-1),dIdt.reshape(-1),dRdt.reshape(-1)))
    
    return dydt


S0 = tot_population.reshape(-1)

I0 = np.zeros((num,num))
I0[int((num-1)/2.),int((num+5)/2.)] = 10
I0 = I0.reshape(-1)

R0 = np.zeros((num,num)).reshape(-1)

y0 = np.concatenate((S0, I0, R0))

"""
print(y0[0:num**2].reshape((num,num)))
print(y0[num**2:2*num**2].reshape((num,num)))
print(y0[2*num**2:3*num**2].reshape((num,num)))
"""

time = np.arange(300)
beta, gamma = 0.02, 0.01

sol = si.odeint(network, y0, time, args=(beta,gamma))
frames =[];
for i in range (len(time)):
    S = sol[i][0:num**2].reshape((num,num))
    I= sol[i][num**2:2*num**2].reshape((num,num))
    R = sol[i][2*num**2:3*num**2].reshape((num,num))
    
   # plt.imshow(I)
   # plt.colorbar()
   # plt.show()
    frames.append(I)


fig = plt.figure()
plt.axis("off")
ax = fig.add_subplot(111)
cv0 = frames[0]
im = ax.imshow(cv0, origin='lower') # Here make an AxesImage rather than contour
cb = fig.colorbar(im)
tx = ax.set_title('Frame 0')

def animate(i):
    arr = frames[i]
    vmax     = np.max(arr)
    vmin     = np.min(arr)
    im.set_data(arr)
    im.set_clim(vmin, vmax)
    tx.set_text('Network Model'.format(i))
    # In this version you don't have to do anything to the colorbar,
    # it updates itself when the mappable it watches (im) changes

ani = animation.FuncAnimation(fig, animate, frames=len(time))

plt.show()

ani.save('network_model_without_fear.mp4')


