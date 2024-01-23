import matplotlib.pyplot as plt 
import numpy as np

"""
2a. Include a hand or computer drawn sketch (doesn't have to be to scale) of the Earth's orbit, the Sun's orbit, and your coordinate axes. 
You can orient the orbit however you prefer.
"""

time = np.arange(0,364,1)

def y_circle_pos(x, radius):
    return np.sqrt(radius**2 - x**2)

def y_circle_neg(x, radius):
    return -np.sqrt(radius**2 - x**2)

e_xorbit = np.linspace(-1, 1, int(len(time)/2))
s_xorbit = np.linspace(-0.1, 0.1, int(len(time)/2))

e_yorbit_pos = y_circle_pos(e_xorbit[::-1], 1)
e_yorbit_neg = y_circle_neg(e_xorbit, 1)
s_yorbit_pos = y_circle_pos(s_xorbit[::-1], 0.1)
s_yorbit_neg = y_circle_neg(s_xorbit, 0.1)

e_yorbit = list(e_yorbit_pos) + list(e_yorbit_neg[1:,])
s_yorbit = list(s_yorbit_pos) + list(s_yorbit_neg[1:,])
e_xorbit = list(reversed(e_xorbit)) + list(e_xorbit[1:,])
s_xorbit = list(reversed(s_xorbit)) + list(s_xorbit[1:,])

plt.scatter([0,1,0.1], [0,0,0], color = 'k')
plt.plot(e_xorbit, e_yorbit, color = 'r')
plt.plot(s_xorbit, s_yorbit, color = 'b')
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("Figures/hw_two_figures/2a.png", bbox_inches = "tight")
plt.show()

"""
2c. Using Euler’s method, we can integrate the positions and velocities:
r’i = ri + vi Δt
v’i = vi + dvi/dt Δt
For this week's assignment, we will set dvi/dt = 0. Using Δt = 1 day, numerically integrate the equations of motion over multiple Earth orbits 
and plot each of the values in the table above as a function of time. Include your plot(s) in your submission. 
"""

def velocity(cur_v, delta_v, delta_time):
    return np.array(cur_v) + delta_v*delta_time

def position(cur_r, cur_v, delta_time):
    return np.array(cur_r) + np.array(cur_v*delta_time)

e_pos = [np.array([1,0])]
e_vel = [np.array([-1,1])]
s_pos = [np.array([0.1,0])]
s_vel = [np.array([-0.1,0.1])]
for i in range(0, len(time)):
    e_pos.append(position(e_pos[i], e_vel[i], 1))
    e_vel.append(velocity(e_vel[i], 0, i))
    s_pos.append(position(s_pos[i], s_vel[i], 1))
    s_vel.append(velocity(s_vel[i], 0, i))

plt.scatter(list(list(zip(*e_pos))[0]), list(list(zip(*e_pos))[1]), color = 'r', label = "Earth")
plt.scatter(list(list(zip(*s_pos))[0]), list(list(zip(*s_pos))[1]), color = 'b', label = "Sun")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.savefig("Figures/hw_two_figures/2c.png", bbox_inches = "tight")
plt.show()

plt.plot(time, list(list(zip(*e_pos))[0])[0:-1], label = "xEarth")
plt.plot(time, list(list(zip(*s_pos))[0])[0:-1], label = "xSun")
plt.plot(time, list(list(zip(*e_pos))[1])[0:-1], label = "yEarth")
plt.plot(time, list(list(zip(*s_pos))[1])[0:-1], label = "ySun")
plt.xlabel("time")
plt.ylabel("equation of motion position variables")
plt.legend()
plt.savefig("Figures/hw_two_figures/2ci.png", bbox_inches = "tight")
plt.show()

plt.plot(time, list(list(zip(*e_vel))[0])[0:-1], label = "Vx, Earth")
plt.plot(time, list(list(zip(*s_vel))[0])[0:-1], label = "Vx, Sun")
plt.plot(time, list(list(zip(*e_vel))[1])[0:-1], label = "Vy, Earth")
plt.plot(time, list(list(zip(*s_vel))[1])[0:-1], label = "Vy, Sun")
plt.xlabel("time")
plt.ylabel("equation of motion velocity variables")
plt.legend()
plt.savefig("Figures/hw_two_figures/2cii.png", bbox_inches = "tight")
plt.show()