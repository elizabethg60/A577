import matplotlib.pyplot as plt
import numpy as np

#info extracted from Table 1 in Mann et al 2013
radius = [0.3863, 0.6398, 0.4840, 0.4183, 0.5477, 0.4712, 0.1869, 0.3924, 0.7949, 0.5773, 0.5673, 0.3982, 0.4546, 0.7390, 0.2990, 0.6697, 0.3561, 0.3232, 0.5472, 0.6611, 0.6010, 0.7784]
radius_err = [0.0021, 0.0046, 0.0084, 0.0070, 0.0048, 0.0086, 0.0012, 0.0033, 0.0062, 0.0131, 0.0137, 0.0091, 0.0182, 0.0190, 0.0100, 0.0089, 0.0039, 0.0061, 0.0067, 0.0048, 0.0072, 0.0053]

mass = [0.405, 0.711, 0.490, 0.403, 0.572, 0.494, 0.159, 0.392, 0.767, 0.630, 0.617, 0.390, 0.447, 0.740, 0.308, 0.749, 0.330, 0.257, 0.573, 0.727, 0.656, 0.771]
mass_err = [0.041, 0.071, 0.049, 0.040, 0.057, 0.049, 0.016, 0.039, 0.124, 0.063, 0.062, 0.039, 0.045, 0.119, 0.031, 0.075, 0.033, 0.026, 0.057, 0.073, 0.066, 0.124]

luminosity = [0.02256, 0.11174, 0.03694, 0.02228, 0.05181, 0.03694, 0.00342, 0.02134, 0.27615, 0.07316, 0.06875, 0.02214, 0.02834, 0.21609, 0.01181, 0.15953, 0.01573, 0.00927, 0.05250, 0.14606, 0.08468, 0.28046]
luminosity_err = [0.00027, 0.00167, 0.00051, 0.00044, 0.00058, 0.00075, 0.00003, 0.00030, 0.00356, 0.00276, 0.00229, 0.00024, 0.00072, 0.00317, 0.00021, 0.00396, 0.00019, 0.00011, 0.00069, 0.00196, 0.00125, 0.00389]

temperature = [3602, 4176, 3646, 3457, 3731, 3695, 3238, 3532, 4704, 3953, 3926, 3537, 3520, 4588, 3487, 4475, 3417, 3142, 3744, 4399, 4025, 4773]
temperature_err = [13, 19, 34, 35, 16, 35, 11, 17, 21, 41, 37, 41, 66, 58, 62, 33, 17, 29, 27, 16, 24, 20]

Fe_H = [-0.3, 0.49, 0.24, -0.31, -0.05, 0.21, -0.06, -0.4, -0.38, -0.28, -0.01, -0.04, -0.37, 0.01, -0.06, -0.15, 0.01, -0.23, -0.3, -0.06, -0.27, -0.22, -0.23]

#prescriptions for radius, mass, and luminosity from Mann et al 2013
def radius_fct(temp):
    return -16.883 + (1.18*10**(-2))*temp - (2.709*10**(-6))*temp**2 + (2.105*10**(-10))*temp**3

def luminosity_fct(temp):
    return -0.781 + (7.4*10**(-4))*temp - (2.49*10**(-7))*temp**2 + (2.95*10**(-11))*temp**3

def mass_fct(temp):
    return -22.297 + (1.544*10**(-2))*temp - (3.488*10**(-6))*temp**2 + (2.65*10**(-10))*temp**3

#prescriptions for radius and luminosity from Boyajian et al 2012
def radius_b12(temp):
    return -10.8828 + (7.18727*10**(-3))*temp - (1.50957*10**(-6))*temp**2 + (1.07572*10**(-10))*temp**3

def luminosity_b12(temp):
    return 10**(-5960.5710 + 4831.6912*np.log10(temp) - 1306.9966*(np.log10(temp))**2 + 117.9716*(np.log10(temp))**3)

#reproduce figure
# # Define a colormap
# cmap = plt.get_cmap('viridis')

# # Normalize the weights to range [0, 1]
# norm = plt.Normalize(min(Fe_H), max(Fe_H))
# normalized_weights = norm(Fe_H)
#c=weights_rescaled, cmap=cmap, 

fig, axs = plt.subplots(3, sharex=True, sharey=False, gridspec_kw={'hspace': 0,})
axs[0].errorbar(temperature, radius, yerr = radius_err, xerr=temperature_err, fmt="o")
axs[0].plot(np.arange(min(temperature), max(temperature)), radius_fct(np.arange(min(temperature), max(temperature))), color = 'b', label = "Our Fit")
axs[0].plot(np.arange(min(temperature), max(temperature)), radius_b12(np.arange(min(temperature), max(temperature))), color = "pink", label = "B12")
axs[0].set_ylabel("Radius ($R_{Sun}$)")
axs[0].legend()
#residuals
axs[1].errorbar(temperature, luminosity, yerr = luminosity_err, xerr=temperature_err, fmt="o")
axs[1].plot(np.arange(min(temperature), max(temperature)), luminosity_fct(np.arange(min(temperature), max(temperature))), color = 'b', label = "Our Fit")
axs[1].plot(np.arange(min(temperature), max(temperature)), luminosity_b12(np.arange(min(temperature), max(temperature))), color = "pink", label = "B12")
axs[1].set_ylabel("Luminosity ($L_{Sun}$)")
#residuals
axs[2].errorbar(temperature, mass, yerr = mass_err, xerr=temperature_err, fmt="o")
axs[2].plot(np.arange(min(temperature), max(temperature)), mass_fct(np.arange(min(temperature), max(temperature))), color = 'b', label = "Our Fit")
axs[2].set_ylabel("Mass ($M_{Sun})$")
axs[2].set_xlabel("$T_{eff}$(K)")
plt.savefig("figure.png")
plt.show()

