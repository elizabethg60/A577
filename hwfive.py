import matplotlib.pyplot as plt 
import numpy as np

mass_sun = 1.99*10**33 #cgs
mass_earth = 5.97*10**27
mass_venus = 4.87*10**27
total_mass = mass_sun + mass_earth + mass_venus
G = 6.67*10**(-8)

def position(cur_r, cur_v, delta_time):
    return np.array(cur_r) + np.array(cur_v*delta_time)

def velocity(cur_v, delta_v, delta_time):
    return np.array(cur_v) + np.array(delta_v*delta_time)

def acceleration(i_r, j_r1, j_r2, mass1, mass2):
    result1 = ((G*mass1*(i_r - j_r1))/((((i_r[0] - j_r1[0])**2) + ((i_r[1] - j_r1[1])**2))**(3/2)))
    result2 = ((G*mass2*(i_r - j_r2))/((((i_r[0] - j_r2[0])**2) + ((i_r[1] - j_r2[1])**2))**(3/2)))

    return -(result1 + result2)

au_cm = 1.496*10**13 #cm in 1 au
venus_cm = 1.077*10**13 #cm in 0.72 au
vel_earth = 2980000 #vel of earth in cm/s
vel_venus = 3500000

step_arr = np.linspace(0.1, 1, 5)
for step in [step_arr[0]]:
    e_pos = [np.array([au_cm,0])] 
    e_vel = [np.array([0, vel_earth])]
    v_pos = [np.array([venus_cm,0])] 
    v_vel = [np.array([0, vel_venus])]
    s_pos = [np.array([-(total_mass - mass_earth*au_cm - mass_venus*venus_cm)/mass_sun,0])]
    s_vel = [np.array([0,-(total_mass - mass_earth*vel_earth - mass_venus*vel_venus)/mass_sun])]

    time = np.arange(0, 10*365*86400, step*86400)
    for i in range(0,len(time)):
        e_pos.append(position(e_pos[i], e_vel[i], step*86400))
        e_vel.append(velocity(e_vel[i], acceleration(e_pos[i],s_pos[i],v_pos[i], mass_sun, mass_venus), step*86400))
        v_pos.append(position(v_pos[i], v_vel[i], step*86400))
        v_vel.append(velocity(v_vel[i], acceleration(v_pos[i],s_pos[i],e_pos[i], mass_sun, mass_earth), step*86400))
        s_pos.append(position(s_pos[i], s_vel[i], step*86400))
        s_vel.append(velocity(s_vel[i], acceleration(s_pos[i],e_pos[i],v_pos[i], mass_earth, mass_venus), step*86400))

    plt.scatter(list(list(zip(*e_pos))[0]), list(list(zip(*e_pos))[1]), color = 'r', label = "Earth", s = 1)
    plt.scatter(list(list(zip(*v_pos))[0]), list(list(zip(*v_pos))[1]), color = 'g', label = "Venus", s = 1)
    plt.scatter(list(list(zip(*s_pos))[0]), list(list(zip(*s_pos))[1]), color = 'b', label = "Sun")
    plt.xlabel("x (cm)")
    plt.ylabel("y (cm)")
    plt.title('Δt = {} days'.format(step))
    plt.legend()
    plt.show()