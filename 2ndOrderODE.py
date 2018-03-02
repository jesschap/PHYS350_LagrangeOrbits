from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


def model(x,t):
  y = x[0]
  dy = x[1]  
  
  xdot = [[],[]]
  xdot[0] = dy
  xdot[1] = np.sin(y)  # equation for y'', replace y' with x[1]
  
  return xdot

conds = [[],[]]
conds[0] = np.pi/9  # initial condition of y
conds[1] = 0        # initial condition of y'
t_f = 15
t_i = 0
time = np.linspace(0,t_f,t_f*100)
z = odeint(model,conds,time)


plt.plot(time[t_i*100:],z[t_i*100:,0],'g:')
#plt.legend(['y (real)', 'y (approx)'])
plt.xlabel('Time')
plt.show()
