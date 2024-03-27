import matplotlib.pyplot as plt
import numpy as np

"""
Create a plot or series of plots to demonstrate the effect of varying the mass of Venus 
-- from 1 lunar mass to 1 Jupiter mass in increments of your choice -- on your O-C plot
"""

#required parameters (all in cgs units)
mass_sun = 1.99*10**33
mass_earth = 5.97*10**27
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

# integrate Sun-Earth-Venus system for 0.1 days with varying Venus mass
mass_venus_values = np.linspace(7.4*10**25, 2*10**30, 5)
trial_values = np.linspace(0.1, 1, 5)
for mass_venus in mass_venus_values:
    total_mass = mass_sun + mass_earth + mass_venus

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
    time = np.arange(0, 51*365*86400, trial_values[0]*86400)
    for i in range(0,len(time)):
        e_pos.append(position(e_pos[i], e_vel[i], trial_values[0]*86400))
        v_pos.append(position(v_pos[i], v_vel[i], trial_values[0]*86400))
        s_pos.append(position(s_pos[i], s_vel[i], trial_values[0]*86400))
        e_vel.append(velocity(e_vel[i], acceleration_3bodies(e_pos[i+1],s_pos[i+1],v_pos[i+1], mass_sun, mass_venus), trial_values[0]*86400))
        v_vel.append(velocity(v_vel[i], acceleration_3bodies(v_pos[i+1],s_pos[i+1],e_pos[i+1], mass_sun, mass_earth), trial_values[0]*86400))
        s_vel.append(velocity(s_vel[i], acceleration_3bodies(s_pos[i+1],e_pos[i+1],v_pos[i+1], mass_earth, mass_venus), trial_values[0]*86400))

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

    transit_number = range(len(transit_time))
    #find line of best fit
    a, b = np.polyfit(transit_number, transit_time, 1)

    #plot transit timing variations using an O-C (observed minus calculated) plot
    observed = transit_time
    calculated = a*transit_number+b
    plt.scatter(observed, observed-calculated, label = f'Venus mass: {mass_venus:.1e} g')

plt.xlabel("O (minutes)")
plt.ylabel("O-C (minutes)")
plt.legend()
plt.savefig("Figures/exam_two_figures/OC_venus_varying.png", bbox_inches = "tight")
plt.show()