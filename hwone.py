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
h = 0.001
for i in t:
    y_euler.append(Euler(t0, y0, h, i, func))

#visuals
plt.plot(t, y, label = "scipy", linewidth = 5)
plt.plot(t, y_euler, label = "euler", color = 'r', linewidth = 5, alpha = 0.5)
plt.xlabel("t", fontsize = 12)
plt.ylabel("z", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.legend(fontsize = 12)
plt.savefig("Figures/hw_one_figures/euler.png", bbox_inches='tight')
plt.show()
