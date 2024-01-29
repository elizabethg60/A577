import matplotlib.pyplot as plt 
import numpy as np

"""
2a. Integrate the Sun-Earth system using Euler’s method for five different trial values of Δt ranging from 0.1 to 1 days. 
Track the relative change in energy and angular momentum.
"""
mass_sun = 2*10**33 #grams
mass_earth = 5.97*10**27
G = (6.67*10**(-8))/((1.496*10**13)**2)

def position(cur_r, cur_v, delta_time):
    return np.array(cur_r) + np.array(cur_v*delta_time)

def velocity(cur_v, delta_v, delta_time):
    return np.array(cur_v) + delta_v*delta_time

def acceleration(i_r, j_r, mass):
    return -np.sum((G*mass*(i_r-j_r))/(((i_r[0] - j_r[0])**2 + (i_r[1] - j_r[1])**2)**(3/2)))

e_pos = [np.array([1,0])]
e_vel = [np.array([-1,1])]
s_pos = [np.array([0.1,0])]
s_vel = [np.array([-0.1,0.1])]

step = 0.1
time = np.arange(0,364,step)
for i in range(0,len(time)):
    e_pos.append(position(e_pos[i], e_vel[i], step))
    e_vel.append(velocity(e_vel[i], acceleration(e_pos[i],s_pos[i], mass_sun), step))
    s_pos.append(position(s_pos[i], s_vel[i], step))
    s_vel.append(velocity(s_vel[i], acceleration(s_pos[i],e_pos[i], mass_earth), step))

"""
2b. For each Δt, plot x vs. y for multiple orbits of the Earth
"""









plt.scatter(list(list(zip(*e_pos))[0]), list(list(zip(*e_pos))[1]), color = 'r', label = "Earth")
plt.scatter(list(list(zip(*s_pos))[0]), list(list(zip(*s_pos))[1]), color = 'b', label = "Sun")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
#plt.savefig("Figures/hw_three_figures/2b.png", bbox_inches = "tight")
plt.show()

# plt.plot(time, list(list(zip(*e_pos))[0])[0:-1], label = "xEarth")
# plt.plot(time, list(list(zip(*s_pos))[0])[0:-1], label = "xSun")
# plt.plot(time, list(list(zip(*e_pos))[1])[0:-1], label = "yEarth")
# plt.plot(time, list(list(zip(*s_pos))[1])[0:-1], label = "ySun")
# plt.xlabel("time")
# plt.ylabel("equation of motion position variables")
# plt.legend()
# plt.savefig("Figures/hw_two_figures/2ci.png", bbox_inches = "tight")
# plt.show()

# plt.plot(time, list(list(zip(*e_vel))[0])[0:-1], label = "Vx, Earth")
# plt.plot(time, list(list(zip(*s_vel))[0])[0:-1], label = "Vx, Sun")
# plt.plot(time, list(list(zip(*e_vel))[1])[0:-1], label = "Vy, Earth")
# plt.plot(time, list(list(zip(*s_vel))[1])[0:-1], label = "Vy, Sun")
# plt.xlabel("time")
# plt.ylabel("equation of motion velocity variables")
# plt.legend()
# plt.savefig("Figures/hw_two_figures/2cii.png", bbox_inches = "tight")
# plt.show()