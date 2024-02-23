import numpy as np
import matplotlib.pyplot as plt

def z(t):
    return ((t-np.pi)**4 + np.sin(t) + 1)

def dzdt(t):
    return 4*(t-np.pi)**3 + np.cos(t)

def dgdt(t):
    return 12*(t-np.pi)**2 - np.sin(t)

def delta_t(t):
    return -dzdt(t)/dgdt(t)

# def D(func, x, h):
# # computes the numerical derivative 
#     return (func(x+h)-func(x))/h

# h = 10**(-6) 
# dxdt = []
# for i in t:
#     dxdt.append(D(z, i, h))

t = np.linspace(0, 6, 100)

guess_time = t[int(len(t)/2)]
guess_dzdt = dzdt(guess_time)
tol = 10**(-8)
while np.abs(0 - guess_dzdt) > tol:
    guess_time = guess_time + delta_t(guess_time)
    guess_dzdt = dzdt(guess_time)


# #determine t nearest to minimum 
# closest_ind = 0
# delta_min = np.abs(0 - dzdt(t[0]))
# for i in range(1, len(t)):
#     if np.abs(0 - dzdt(t[i])) < delta_min:
#         delta_min = np.abs(0 - dzdt(t[i]))
#         closest_ind = i

# tol = 10**(-16)
# guess = t[0]
# min_guess = dzdt(guess)
# while np.abs(0 - guess) > tol: 
#     guess = guess + delta_t()


plt.plot(t, z(t), label = 'z')
# plt.scatter(t, dxdt, label = 'dz')
plt.scatter(guess_time, guess_dzdt)
# plt.scatter(t[closest_ind], dzdt(t[closest_ind]))
#plt.scatter(t, dzdt(t), label = 'dz/dt')
# plt.scatter(t, dgdt(t), label = 'dg/dt')
plt.legend()
plt.show()

print(guess_time, guess_dzdt)