import matplotlib.pyplot as plt 
from A577_package import Euler
import numpy as np

#define function
def func(y, t):
    return np.sqrt(y)
y0 = 1            
t0 = 0              
tf = 5
t = np.arange(t0, tf, 0.01)  

#initially solving with numerical integration package to compare answers
from scipy.integrate import odeint
y = odeint(func, y0, t) 

#my code for solving ODE
y_euler = []
h = 0.00001
for i in t:
    y_euler.append(Euler(t0, y0, h, i, func))

#visuals
plt.plot(t, y, label = "scipy")
plt.plot(t, y_euler, label = "euler")
plt.legend()
plt.show()