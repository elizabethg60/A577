import matplotlib.pyplot as plt 
import numpy as np

"""
2a. Integrate the Sun-Earth system using Euler’s method for five different trial values of Δt ranging from 0.1 to 1 days. 
Track the relative change in energy and angular momentum.

2b. For each Δt, plot x vs. y for multiple orbits of the Earth
"""

mass_sun = 1.99*10**33 #cgs
mass_earth = 5.97*10**27
radius_sun = 6.96*10**10
radius_earth = 6.37*10**8
ang_v_earth =(2*np.pi)/(1*24*3600)
ang_v_sun = (2*np.pi)/(25*24*3600)
G = 6.67*10**(-8)

def position(cur_r, cur_v, delta_time):
    return np.array(cur_r) + np.array(cur_v*delta_time)

def velocity(cur_v, delta_v, delta_time):
    return np.array(cur_v) + np.array(delta_v*delta_time)

def acceleration(i_r, j_r, mass):
    return -((G*mass*(i_r - j_r))/((((i_r[0] - j_r[0])**2) + ((i_r[1] - j_r[1])**2))**(3/2)))

def energy(e_pos, e_vel, s_pos, s_vel):
    return -((G*mass_sun*mass_earth)/(np.abs(np.linalg.norm(e_pos-s_pos)))) \
    + 0.5*(mass_earth)*((np.linalg.norm(e_vel))**2) + 0.5*(mass_sun)*((np.linalg.norm(s_vel))**2)

def momentum(mass, cur_r, cur_v):
    return mass*np.linalg.norm(cur_r)*np.linalg.norm(cur_v)

au_cm = 1.496*10**13 #cm in 1 au
vel = 2980000 #vel of earth in cm/s

step_arr = np.linspace(0.1, 1, 5)
for step in step_arr:
    e_pos = [np.array([au_cm,0])] 
    e_vel = [np.array([0, vel])]
    s_pos = [np.array([-(mass_earth*au_cm)/mass_sun,0])]
    s_vel = [np.array([0,-(mass_earth*vel)/mass_sun])]

    time = np.arange(0, 10*365*86400, step*86400)
    for i in range(0,len(time)):
        e_pos.append(position(e_pos[i], e_vel[i], step*86400))
        e_vel.append(velocity(e_vel[i], acceleration(e_pos[i],s_pos[i], mass_sun), step*86400))
        s_pos.append(position(s_pos[i], s_vel[i], step*86400))
        s_vel.append(velocity(s_vel[i], acceleration(s_pos[i],e_pos[i], mass_earth), step*86400))

    plt.scatter(list(list(zip(*e_pos))[0]), list(list(zip(*e_pos))[1]), color = 'r', label = "Earth", s = 1)
    plt.scatter(list(list(zip(*s_pos))[0]), list(list(zip(*s_pos))[1]), color = 'b', label = "Sun")
    plt.xlabel("x (cm)")
    plt.ylabel("y (cm)")
    plt.title('Δt = {} days'.format(step))
    plt.legend()
    plt.savefig("Figures/hw_three_figures/2_step_{}.png".format(step), bbox_inches = "tight")
    plt.show()

step_arr = np.linspace(0.1, 1, 5)
for step in step_arr:
        e_pos = [np.array([au_cm,0])] 
        e_vel = [np.array([0, vel])]
        s_pos = [np.array([-(mass_earth*au_cm)/mass_sun,0])]
        s_vel = [np.array([0,-(mass_earth*vel)/mass_sun])]

        energy_arr = []
        e_momentum = []
        for i in range(0,len(time)):
            e_pos.append(position(e_pos[i], e_vel[i], step*86400))
            e_vel.append(velocity(e_vel[i], acceleration(e_pos[i],s_pos[i], mass_sun), step*86400))
            s_pos.append(position(s_pos[i], s_vel[i], step*86400))
            s_vel.append(velocity(s_vel[i], acceleration(s_pos[i],e_pos[i], mass_earth), step*86400))
            energy_arr.append(energy(e_pos[i], e_vel[i], s_pos[i], s_vel[i]))

        diff_list_energy = []
        for x in range(0,len(energy_arr)):
            diff_list_energy.append((energy_arr[x] - energy_arr[0])/energy_arr[0])

        plt.scatter((time/86400), diff_list_energy, label = "{} delta time".format(step))
plt.xlabel("time (days)")
plt.ylabel("energy")
plt.title("relative change in energy")
plt.legend()
plt.savefig("Figures/hw_three_figures/2_energy.png", bbox_inches = "tight")
plt.show()

for step in step_arr:
        e_pos = [np.array([au_cm,0])] 
        e_vel = [np.array([0, vel])]
        s_pos = [np.array([-(mass_earth*au_cm)/mass_sun,0])]
        s_vel = [np.array([0,-(mass_earth*vel)/mass_sun])]

        energy_arr = []
        e_momentum = []
        for i in range(0,len(time)):
            e_pos.append(position(e_pos[i], e_vel[i], step*86400))
            e_vel.append(velocity(e_vel[i], acceleration(e_pos[i],s_pos[i], mass_sun), step*86400))
            s_pos.append(position(s_pos[i], s_vel[i], step*86400))
            s_vel.append(velocity(s_vel[i], acceleration(s_pos[i],e_pos[i], mass_earth), step*86400))
            e_momentum.append(momentum(mass_earth, e_pos[i], e_vel[i]) + momentum(mass_sun, s_pos[i], s_vel[i]))

        diff_list_mom = []
        for x in range(0,len(e_momentum)):
            diff_list_mom.append((e_momentum[x] - e_momentum[0])/e_momentum[0])

        plt.scatter((time/86400), diff_list_mom, label = "{} delta time".format(step))
plt.xlabel("time (days)")
plt.ylabel("angular momentum")
plt.title('relative change in angular momentum')
plt.legend()
plt.savefig("Figures/hw_three_figures/2_momentum", bbox_inches = "tight")
plt.show()