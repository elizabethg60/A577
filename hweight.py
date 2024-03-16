import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

"""
We plot transit timing variations using an O-C (observed minus calculated) plot, 
where calculated refers to the times from a best-fit linear ephemeris. 
In our case, the "observed" transit times are those computed by your simulation code. 
You can fit a line to those transit times to get a best-fit linear ephemeris, 
which I recommend parametrizing transit time as a function of transit number. 
"""

#required parameters (all in cgs units)
mass_sun = 1.99*10**33
mass_earth = 5.97*10**27
mass_venus = 4.87*10**27 #0.0
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
        if e_pos[i][1] < 0 and e_pos[i+1][1] > 0:
            e_pos_y = e_pos[i]
            e_vel_y = e_vel[i] 
            v_pos_y = v_pos[i]
            v_vel_y = v_vel[i]
            s_pos_y = s_pos[i]
            s_vel_y = s_vel[i]

            #We want to minimize |rs|, which occurs when: g = rs •drs/dt
            #We can refine our mid transit time solution using: ∆ t = - g (∂g/∂t)-1

            min_time = time[i]
            min_value = e_pos_y[1] * e_vel_y[1]
            ay = acceleration_3bodies(e_pos[i+1],s_pos[i+1],v_pos[i+1], mass_sun, mass_venus)[1]
            dg_dt = e_vel_y[1]**2 + e_pos_y[1] * ay

            tol = 10**(-4)
            while np.abs(min_value - 0) > tol:
                delta_t = -min_value/dg_dt
                min_time = min_time + delta_t

                e_pos_y = position(e_pos_y, e_vel_y, delta_t)
                v_pos_y = position(v_pos_y, v_vel_y, delta_t)
                s_pos_y = position(s_pos_y, s_vel_y, delta_t)
                e_vel_y = velocity(e_vel_y, acceleration_3bodies(e_pos_y,s_pos_y,v_pos_y, mass_sun, mass_venus), delta_t)
                v_vel_y = velocity(v_vel_y, acceleration_3bodies(v_pos_y,s_pos_y,e_pos_y, mass_sun, mass_earth), delta_t)
                s_vel_y = velocity(s_vel_y, acceleration_3bodies(s_pos_y,e_pos_y,v_pos_y, mass_earth, mass_venus), delta_t)

                min_value = e_pos_y[1] * e_vel_y[1]
                ay = acceleration_3bodies(e_pos_y,s_pos_y,v_pos_y, mass_sun, mass_venus)[1]
                dg_dt = e_vel_y[1]**2 + e_pos_y[1] * ay
            
            transit_pos.append(e_pos_y)
            transit_time.append(min_time/(60))
            transit_sep.append(e_pos_y[1]/(1.496*10**13))

    #plot transit number vs transit time
    transit_number = range(len(transit_time))
    #find line of best fit
    a, b = np.polyfit(transit_number, transit_time, 1)
    plt.scatter(transit_number, transit_time, color = 'k')
    plt.plot(transit_number, a*transit_number+b)
    plt.xlabel("transit number")
    plt.ylabel("transit time (minutes)")
    plt.title('Δt = {} days'.format(step))
    plt.savefig("Figures/hw_eight_figures/transit_linear_venus.png", bbox_inches = "tight")
    plt.show()

    #plot transit timing variations using an O-C (observed minus calculated) plot
    observed = transit_time
    calculated = a*transit_number+b

    """
    Create an O-C plot from your simulation with Venus' mass set to zero. 
    The x-axis should be the time of mid-transit from your simulation, and the y-axis should be  
    mid-transit from your simulation minus the best fit linear ephemeris. Plot the y-axis values in minutes. 
    The values on the y-axis should be quite small (<1 min).

    Create an O-C plot from your simulation with Venus' mass set to its true value 
    The values on the y-axis should be larger but still small (a few minutes).
    """
    plt.scatter(observed, observed-calculated, color = 'k')
    plt.xlabel("O (minutes)")
    plt.ylabel("O-C (minutes)")
    plt.savefig("Figures/hw_eight_figures/OC_venus.png", bbox_inches = "tight")
    plt.show()