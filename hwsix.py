import matplotlib.pyplot as plt
import numpy as np

"""
2a) We can find the minimum by finding the t for which : g = dz/dt = 0
Create an array of t and plot z vs. t.
Hone in on the t nearest to the minimum. Refine your solution using ∆ t = - g (dg/dt)-1
with multiple iterations until you reach a precise solution. 
Add your more precise value for the t of the minimum as a labeled point on your plot and state it to several decimal places.
"""

def z(t):
# returns z value for a given t
    return (t-np.pi)**4 + np.sin(t) + 1

def g(t):
# returns g value for a given t
    return 4*(t-np.pi)**3 + np.cos(t)

def dgdt(t):
# returns dgdt for a given t
    return 12*(t-np.pi)**2 - np.sin(t)

def delta_t(t):
# returns delta t for a given t
    return -g(t)/dgdt(t)

#array of t values
t = np.linspace(0, 6, 100)

#iterate using delta t to find min of z func
#initial guess
min_t = t[int(len(t)/2)]
min_value = g(min_t)

tol = 10**(-4)
while np.abs(min_value - 0) > tol:
    min_t += delta_t(min_t)
    min_value = g(min_t)

#plot
plt.plot(t, z(t), label = 'z')
plt.scatter(min_t, min_value)
plt.xlabel("t")
plt.ylabel("z")
plt.legend()
plt.savefig("Figures/hw_six_figures/min_z.png", bbox_inches = "tight")
plt.show()

print("min z occurs at t: {}".format(min_t))
print("min z: {}".format(min_value))