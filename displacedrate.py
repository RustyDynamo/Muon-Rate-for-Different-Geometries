# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import quad
!pip install daemonflux
import daemonflux
from daemonflux import Flux

"""# Rate through two displaced areas

In this section, we look at the muon rate passing through two circular areas that are both horizontally and vertically displaced. This is meant to simulate the rate of muons passing through two modules situated on adjacent strings and situated one below the other. Here, the two areas are flat circles whose area is the circular projection of an mDOM.
When plotting, we can vary the vertical and the horizontal displacement by keeping the other constant. Naturally, both plots should have the same muon rate at a same point (e.g. the real separation of modules in the Upgrade).
"""

# Constants
r = 17.8                   # cm
area = np.pi * r**2        # cm^2
# Values of mDOM spacing in ice for IceCube Upgrade
default_x = 4500            # cm
default_z = 300            # cm

# momentum grid
p_min = 10.0               # GeV
pgrid = np.linspace(p_min, 10000, 1000)
dp    = pgrid[1] - pgrid[0] #grid steps

# number of angle grid steps
n_theta = 200

# initializing daemonflux
flux_model = Flux(location='generic')

#We first normalize the muon flux and integrate over the momentum grid
def S_of_theta(theta_rad):
    """
    S(θ) = ∫[Φ(p,θ) / p^3] dp, from p_min → 1e4 GeV
    """
    theta_deg = np.degrees(theta_rad) #conversion of angle from radians to degrees
    Φ = flux_model.flux(pgrid, theta_deg, 'muflux') #muon flux function
    return np.sum(Φ / (pgrid**3)) * dp

#Next, we calculate the rate for two disks sparated by specific horizontal and vertical separations (x,z).
#This is done by integrating over all possible angles that still hit both areas.
def rate_quad(x, z):
    H = np.hypot(x, z) #Hypotenuse between centers of areas
    theta_max = np.arctan(r / H) #maximum angle that hits both areas
    thetas = np.linspace(0, theta_max, n_theta)
    dtheta = thetas[1] - thetas[0]

    total = 0.0 #initializing theta integration
    for θ in thetas: #Riemann Integration over solid angle
        Sθ     = S_of_theta(θ)
        Aproj  = area * np.cos(θ) #Area depending on muon angle
        dΩ     = 2*np.pi * np.sin(θ) * dtheta #solid angle element
        total += Sθ * Aproj * dΩ

    return total

#Reference rate at IceCube Upgrade dimensions
ref_rate = rate_quad(default_x, default_z)
print(f"Reference rate (p>10 GeV): {ref_rate:.2e} 1/s")

#Plot of rate as a function of vertical displacement z
z_vals = np.linspace(50, 1000, 40) #values of vertical displacement
z_rates = [rate_quad(default_x, z) for z in z_vals] #associated muon rate for each z value with constant x

plt.figure(figsize=(8,6))
plt.plot(z_vals, z_rates, 'b-o', label='Muon Rate')
plt.axvline(default_z, color='r', linestyle='--', label=f'z = {default_z} cm')
plt.axhline(ref_rate, color='g', linestyle='--', label=f'Rate = {ref_rate:.2e} 1/s')
plt.xlabel("Vertical displacement z (cm)")
plt.ylabel("Muon rate (1/s)")
plt.title("mDOM Muon rate vs. vertical displacement (x = 4500 cm)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

#Plot of rate as a function of horizontal displacement
x_vals = np.linspace(50, 5000, 50) #values of horizontal displacement
x_rates = [rate_quad(x, default_z) for x in x_vals] #associated muon rate for each x value with constant z

plt.figure(figsize=(8,6))
plt.plot(x_vals, x_rates, 'm-o', label='Muon Rate')
plt.axvline(default_x, color='r', linestyle='--', label=f'x = {default_x} cm')
plt.axhline(ref_rate, color='g', linestyle='--', label=f'Rate = {ref_rate:.2e} 1/s')
plt.xlabel("Horizontal displacement x (cm)")
plt.ylabel("Muon rate (1/s)")
plt.title("mDOM Muon rate vs. horizontal displacement (z = 300 cm)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

"""**Note:** The time for the plots might be rather long as the code calculates the muon rate for each individual chosen separation individually. Plotting more points thus drastically increases computation time"""
