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

def acceleration_3bodies(i_r, j_r1, j_r2, mass1, mass2):
# returns current acceleration provided bodies in question position and mass
    result1 = ((G*mass1*(i_r - j_r1))/((((i_r[0] - j_r1[0])**2) + ((i_r[1] - j_r1[1])**2))**(3/2)))
    result2 = ((G*mass2*(i_r - j_r2))/((((i_r[0] - j_r2[0])**2) + ((i_r[1] - j_r2[1])**2))**(3/2)))

    return -(result1 + result2)

def g(y, v_y):
# returns g value 
    return y*v_y

def dgdt(y, v_y, a_y):
# returns dgdt 
    return v_y**2 + y*a_y

def delta_t(y, v_y, a_y):
# returns delta t 
    return -g(y, v_y)/dgdt(y, v_y, a_y)

"""
Using your orbital integrator with a system containing the Sun, the Earth, and Venus. 
At the approximate times of Earth's transit over 50 orbits, compute the magnitude of r_s, 
the projected separation between the center of the Sun and the center of the Earth. Create a plot of your results 
(r_s vs. time, 50 points each at the approximate time of transit). Make sure to label your units.
You'll need to consider when the Earth is transiting in your coordinate system. 
"""

# integrate Sun-Earth-Venus system for 0.1 days
trial_values = np.linspace(0.1, 1, 5)
for step in [trial_values[0]]:
    # hardwire initial positions according to Earth's and Venus's semi-major axis and speed 
        #so Sun's values gives center of mass of 0.0
    e_pos = [np.array([earth_a,0])] 
    e_vel = [np.array([0, earth_vel])]
    v_pos = [np.array([venus_a,0])] 
    v_vel = [np.array([0, venus_vel])]
    s_pos = [np.array([-(total_mass - mass_earth*earth_a - mass_venus*venus_a)/mass_sun,0])]
    s_vel = [np.array([0,-(total_mass - mass_earth*earth_vel - mass_venus*venus_vel)/mass_sun])]

    #collect neccessary values for transit 
    transit_pos = []
    transit_time = []
    transit_sep = []
    #iterate through 50 orbits for given delta time
    time = np.arange(0, 51*365*86400, step*86400)
    for i in range(0,len(time)):
        e_pos.append(position(e_pos[i], e_vel[i], step*86400))
        v_pos.append(position(v_pos[i], v_vel[i], step*86400))
        s_pos.append(position(s_pos[i], s_vel[i], step*86400))
        e_vel.append(velocity(e_vel[i], acceleration_3bodies(e_pos[i+1],s_pos[i+1],v_pos[i+1], mass_sun, mass_venus), step*86400))
        v_vel.append(velocity(v_vel[i], acceleration_3bodies(v_pos[i+1],s_pos[i+1],e_pos[i+1], mass_sun, mass_earth), step*86400))
        s_vel.append(velocity(s_vel[i], acceleration_3bodies(s_pos[i+1],e_pos[i+1],v_pos[i+1], mass_earth, mass_venus), step*86400))

        #look for transit with observer along position x-axis ie change in sign of y value from neg to pos
        min_t = time[int(len(time)/2)]
        min_value = g(min_t)

        tol = 10**(-4)
        while np.abs(min_value - 0) > tol:
            min_t += delta_t(min_t)
            min_value = g(min_t)


        # if e_pos[i][1] < 0 and e_pos[i+1][1] > 0:
        #     r_first = np.abs(e_pos[i][1] - s_pos[i][1])
        #     r_second = np.abs(e_pos[i+1][1] - s_pos[i+1][1])
        #     if r_first > r_second:
        #         transit_pos.append(e_pos[i+1])
        #         transit_time.append(time[i+1]/(3.154*10**7))
        #         transit_sep.append(r_second/(1.496*10**13))
        #     if r_first < r_second:
        #         transit_pos.append(e_pos[i])
        #         transit_time.append(time[i]/(3.154*10**7))
        #         transit_sep.append(r_first/(1.496*10**13))

    # #plot x vs y for given delta time
    # plt.scatter(list(list(zip(*e_pos))[0]), list(list(zip(*e_pos))[1]), color = 'r', label = "Earth")
    # plt.scatter(list(list(zip(*transit_pos))[0]), list(list(zip(*transit_pos))[1]), color = 'g', label = "Earth Transit")
    # plt.scatter(list(list(zip(*v_pos))[0]), list(list(zip(*v_pos))[1]), color = 'b', label = "Venus")
    # plt.scatter(list(list(zip(*s_pos))[0]), list(list(zip(*s_pos))[1]), color = 'k', label = "Sun")
    # plt.xlabel("x (cm)")
    # plt.ylabel("y (cm)")
    # plt.title('Δt = {} days'.format(step))
    # plt.legend()
    # #plt.savefig("Figures/hw_five_figures/transit_xy_{}.png".format(step), bbox_inches = "tight")
    # plt.show()

    #plot r_s vs time
    plt.scatter(transit_time, transit_sep, color = 'k')
    plt.xlabel("time (years)")
    plt.ylabel("Sun to Earth projected separation (AU)")
    plt.title('Δt = {} days'.format(step))
    #plt.savefig("Figures/hw_five_figures/transit_sep_{}.png".format(step), bbox_inches = "tight")
    plt.show()