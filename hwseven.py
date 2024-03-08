import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

"""
Use the algorithm above to hone in on the values of t that minimize rs 
and create an updated plot of rs vs. time of transit for a system with Venus' mass set to 0. 
Your rs values should be smaller than on Homework 5.
"""

#required parameters (all in cgs units)
mass_sun = 1.99*10**33
mass_earth = 5.97*10**27
mass_venus = 0.0
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
            transit_time.append(min_time/(3.154*10**7))
            transit_sep.append(e_pos_y[1]/(1.496*10**13))

    #plot x vs y for given delta time
    plt.scatter(list(list(zip(*e_pos))[0]), list(list(zip(*e_pos))[1]), color = 'r', label = "Earth")
    plt.scatter(list(list(zip(*transit_pos))[0]), list(list(zip(*transit_pos))[1]), color = 'g', label = "Earth Transit")
    plt.scatter(list(list(zip(*v_pos))[0]), list(list(zip(*v_pos))[1]), color = 'b', label = "Venus")
    plt.scatter(list(list(zip(*s_pos))[0]), list(list(zip(*s_pos))[1]), color = 'k', label = "Sun")
    plt.xlabel("x (cm)")
    plt.ylabel("y (cm)")
    plt.title('Δt = {} days'.format(step))
    plt.legend()
    plt.savefig("Figures/hw_seven_figures/transit_xy_{}.png".format(step), bbox_inches = "tight")
    plt.show()

    #plot r_s vs time
    plt.scatter(transit_time, transit_sep, color = 'k')
    plt.xlabel("time (years)")
    plt.ylabel("Sun to Earth projected separation (AU)")
    plt.title('Δt = {} days'.format(step))
    plt.yscale("log")
    plt.savefig("Figures/hw_seven_figures/transit_sep_{}.png".format(step), bbox_inches = "tight")
    plt.show()

"""
Calculate the spacing in mutual Hill radii for each pair of adjacent planets in the TRAPPIST-1 system.
"""

#mass in earth mass
star_mass = 29639.693
b_mass = 1.02 
c_mass = 1.16
d_mass = 0.297
e_mass = 0.772
f_mass = 0.934 
g_mass = 1.148
h_mass = 0.331

#semi-major axis in AU
b_axis = 0.01150
c_axis = 0.01576
d_axis = 0.02219
e_axis = 0.02916
f_axis = 0.03836
g_axis = 0.0467 
h_axis = 0.0617

planets = ["b","c","d","e","f","g","h"]
planets_mass = [b_mass, c_mass, d_mass, e_mass, f_mass, g_mass, h_mass]
planets_axis = [b_axis, c_axis, d_axis, e_axis, f_axis, g_axis, h_axis]

df = pd.DataFrame(columns = ['Planet 1', 'Planet 2', 'mutual Hill radii', 'delta'])

def mutual_radii(a1, a2, m1, m2):
    return ((a1+a2)/2) * ((m1+m2)/(3*star_mass))**(1/3)

def delta(a1, a2, radii):
    return (a2 - a1)/radii

for i in range(0, len(planets_mass)-1):
    current_radii = mutual_radii(planets_axis[i], planets_axis[i+1], planets_mass[i], planets_mass[i+1])
    current_delta = delta(planets_axis[i], planets_axis[i+1], current_radii)
    df = df._append({'Planet 1' : planets[i], 'Planet 2': planets[i+1], 'mutual Hill radii' : current_radii, 'delta' : current_delta},ignore_index = True)

print(df)

"""
Create a cartoon plot of C / O vs. C / H, labeling the regions where we expect to find planets with 
different histories. Each region should be labeled with a short phrase to describe the formation history.
"""
#points extracted using web app from figure 3 in paper 
c_o = [1.97, 1.62, 1.32, 1.097, 0.997, 0.879, 0.768, 0.688]
c_h = [0.386, 0.529, 0.722, 0.894, 0.998, 1.14, 1.28, 1.426]
#create cartoon plot
plt.scatter(c_o, c_h, color = 'k')
plt.plot(c_o, c_h, color = 'k')
plt.axvline(x=1.2)
plt.axvline(x=0.9)
plt.text(1.5, 1.0, 'Super-Stellar:')
plt.text(1.3, 0.9, 'gas accretion close to CO or CO2 ')
plt.text(1.25, 0.8, 'snowlines and outside water snowline')
plt.text(0.95, 1.4, 'Stellar:')
plt.text(0.9, 1.3, 'gravitational collapse or')
plt.text(0.9, 1.2, 'gas accretion followed')
plt.text(0.9, 1.1, 'by inward migration')
plt.text(0.65, 0.8, 'Sub-Stellar:')
plt.text(0.65, 0.7, 'gravitational collapse')
plt.text(0.65, 0.6, 'or core accretion ')
plt.xlabel("C/O")
plt.ylabel("C/H")
plt.savefig("Figures/hw_seven_figures/composition.png", bbox_inches = "tight")
plt.show()
