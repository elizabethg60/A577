import matplotlib.pyplot as plt 
import numpy as np

#required parameters (all in cgs units)
mass_sun = 1.99*10**33
mass_earth = 5.97*10**27
mass_venus = 4.87*10**27
total_mass = mass_sun + mass_earth + mass_venus
G = 6.67*10**(-8)
earth_a = 1.496*10**13 #semi-major axis
venus_a = 1.077*10**13
earth_vel = 2980000
venus_vel = 3500000

def position(cur_r, cur_v, delta_time):
# returns following position provided current position / velocity and delta time 
    return np.array(cur_r) + np.array(cur_v*delta_time)

def velocity(cur_v, delta_v, delta_time):
# returns following velocity provided current velocity, acceleration and delta time 
    return np.array(cur_v) + np.array(delta_v*delta_time)

def acceleration(i_r, j_r, mass):
# returns current acceleration provided bodies in question position and mass
    return -((G*mass*(i_r - j_r))/((((i_r[0] - j_r[0])**2) + ((i_r[1] - j_r[1])**2))**(3/2)))

def energy(e_pos, e_vel, s_pos, s_vel):
# returns total energy (PE + KE of each body) 
    return -((G*mass_sun*mass_earth)/(np.abs(np.linalg.norm(e_pos-s_pos)))) \
    + 0.5*(mass_earth)*((np.linalg.norm(e_vel))**2) + 0.5*(mass_sun)*((np.linalg.norm(s_vel))**2)

def momentum(mass, cur_r, cur_v):
# returns momentum for a given object
    return mass*np.linalg.norm(cur_r)*np.linalg.norm(cur_v)

"""
Integrate the Sun-Earth system using the symplectic method for five different trial values of Δt ranging from 0.1 to 1 days.
Track the fractional relative change in energy and angular momentum. 
Your submission should include ten plots of relative change vs. time in total (or fewer if you combine them into fewer 
plots). 
For each Δt, plot x vs. y for multiple orbits of the Earth.
"""

#integrate Sun-Earth system for 5 different delta time from 0.1 to 1 days
trial_values = np.linspace(0.1, 1, 5)
#track fractional relative change in energy and angular momentum
final_energy_change = []
final_momentum_change = []

for step in trial_values:
    # hardwire initial positions according to Earth's semi-major axis and speed so Sun's values give center of mass of 0.0
    e_pos = [np.array([earth_a,0])] 
    e_vel = [np.array([0, earth_vel])]
    s_pos = [np.array([-(mass_earth*earth_a)/mass_sun,0])]
    s_vel = [np.array([0,-(mass_earth*earth_vel)/mass_sun])]

    #iterate through 10 orbits for given delta time
    energy_arr = []
    momentum_arr = []
    time = np.arange(0, 10*365*86400, step*86400)
    for i in range(0,len(time)):
        e_pos.append(position(e_pos[i], e_vel[i], step*86400))
        s_pos.append(position(s_pos[i], s_vel[i], step*86400))
        e_vel.append(velocity(e_vel[i], acceleration(e_pos[i+1],s_pos[i+1], mass_sun), step*86400))
        s_vel.append(velocity(s_vel[i], acceleration(s_pos[i+1],e_pos[i+1], mass_earth), step*86400))
        energy_arr.append(energy(e_pos[i], e_vel[i], s_pos[i], s_vel[i]))
        momentum_arr.append(momentum(mass_earth, e_pos[i], e_vel[i]) + momentum(mass_sun, s_pos[i], s_vel[i]))

    #calculate relative change in energy and momentum
    energy_change = []
    for x in range(0,len(energy_arr)):
        energy_change.append((energy_arr[x] - energy_arr[0])/energy_arr[0])
    final_energy_change.append(energy_change)

    momentum_change = []
    for x in range(0,len(momentum_arr)):
        momentum_change.append((momentum_arr[x] - momentum_arr[0])/momentum_arr[0])
    final_momentum_change.append(momentum_change)

    #plot x vs y for given delta time
    plt.scatter(list(list(zip(*e_pos))[0]), list(list(zip(*e_pos))[1]), color = 'r', label = "Earth")
    plt.scatter(list(list(zip(*s_pos))[0]), list(list(zip(*s_pos))[1]), color = 'k', label = "Sun")
    plt.xlabel("x (cm)")
    plt.ylabel("y (cm)")
    plt.title('Δt = {} days'.format(step))
    plt.legend()
    plt.savefig("Figures/exam_one_figures/2_step_{}.png".format(step), bbox_inches = "tight")
    plt.show()

#plot relative change in energy and momentum vs time 
for i in range(0,len(final_energy_change)):
    time = np.arange(0, 10*365*86400, trial_values[i]*86400)
    plt.plot((time/86400), np.abs(final_energy_change[i]), label = "Δt = {} days".format(trial_values[i]))
plt.xlabel("time (days)")
plt.ylabel("relative change in energy")
plt.legend()
plt.savefig("Figures/exam_one_figures/2_energy.png", bbox_inches = "tight")
plt.show()


for i in range(0,len(final_momentum_change)):
    time = np.arange(0, 10*365*86400, trial_values[i]*86400)
    plt.plot((time/86400), final_momentum_change[i], label = "Δt = {} days".format(trial_values[i]))
plt.xlabel("time (days)")
plt.ylabel("relative change in angular momentum'")
plt.legend()
plt.savefig("Figures/exam_one_figures/2_momentum", bbox_inches = "tight")
plt.show()

"""
Add Venus to the system. Then plot x vs. y for multiple orbits of the Earth and Venus, Δt = 0.1 day. 
"""

def acceleration_3bodies(i_r, j_r1, j_r2, mass1, mass2):
# returns current acceleration provided bodies in question position and mass
    result1 = ((G*mass1*(i_r - j_r1))/((((i_r[0] - j_r1[0])**2) + ((i_r[1] - j_r1[1])**2))**(3/2)))
    result2 = ((G*mass2*(i_r - j_r2))/((((i_r[0] - j_r2[0])**2) + ((i_r[1] - j_r2[1])**2))**(3/2)))

    return -(result1 + result2)

# integrate Sun-Earth-Venus system for 0.1 days
for step in [trial_values[0]]:
    # hardwire initial positions according to Earth's and Venus's semi-major axis and speed 
        #so Sun's values gives center of mass of 0.0
    e_pos = [np.array([earth_a,0])] 
    e_vel = [np.array([0, earth_vel])]
    v_pos = [np.array([venus_a,0])] 
    v_vel = [np.array([0, venus_vel])]
    s_pos = [np.array([-(total_mass - mass_earth*earth_a - mass_venus*venus_a)/mass_sun,0])]
    s_vel = [np.array([0,-(total_mass - mass_earth*earth_vel - mass_venus*venus_vel)/mass_sun])]

    #iterate through 10 orbits for given delta time
    time = np.arange(0, 10*365*86400, step*86400)
    for i in range(0,len(time)):
        e_pos.append(position(e_pos[i], e_vel[i], step*86400))
        v_pos.append(position(v_pos[i], v_vel[i], step*86400))
        s_pos.append(position(s_pos[i], s_vel[i], step*86400))
        e_vel.append(velocity(e_vel[i], acceleration_3bodies(e_pos[i+1],s_pos[i+1],v_pos[i+1], mass_sun, mass_venus), step*86400))
        v_vel.append(velocity(v_vel[i], acceleration_3bodies(v_pos[i+1],s_pos[i+1],e_pos[i+1], mass_sun, mass_earth), step*86400))
        s_vel.append(velocity(s_vel[i], acceleration_3bodies(s_pos[i+1],e_pos[i+1],v_pos[i+1], mass_earth, mass_venus), step*86400))

    #plot x vs y for given delta time
    plt.scatter(list(list(zip(*e_pos))[0]), list(list(zip(*e_pos))[1]), color = 'r', label = "Earth")
    plt.scatter(list(list(zip(*v_pos))[0]), list(list(zip(*v_pos))[1]), color = 'b', label = "Venus")
    plt.scatter(list(list(zip(*s_pos))[0]), list(list(zip(*s_pos))[1]), color = 'k', label = "Sun")
    plt.xlabel("x (cm)")
    plt.ylabel("y (cm)")
    plt.title('Δt = {} days'.format(step))
    plt.legend()
    plt.savefig("Figures/exam_one_figures/venus_orbit_{}.png".format(step), bbox_inches = "tight")
    plt.show()
    