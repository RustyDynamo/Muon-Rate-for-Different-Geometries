# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import quad
!pip install daemonflux
import daemonflux
from daemonflux import Flux

"""#Rate through two superposed areas

To determine that a muon is truly passing through an mDOM, and not some random event that would trigger one PMT such as radiation, we consider that at least two PMTs would pick up its signal. Thus we calculate the rate of muons passing through 2 areas, starting here with them being directly superposed.
The constants used in this code represent the circular projection of the Hamamatsu R15458-02 Photomultiplier tubes (PMTs) used in the mDOMs of the IceCube Upgrade
"""

# Constants (using dimensions of 2 diametrically opposite PMTs within an mDOM)
diameter = 8.0  # cm
radius = 4.0  # cm
area_circle = np.pi * radius ** 2  # cm^2
distance_between_areas = 41.0  # cm

# Initialize daemonflux
daemonflux = Flux(location='generic')

# Muon flux function with theta dependency
def muon_flux(p, theta):
    return daemonflux.flux(p, theta * 180 / np.pi, 'muflux')  # Convert theta to degrees for DAEMONFLUX

# Effective area function
def effective_area(theta):
    return area_circle * np.cos(theta)  # Projection of area due to angle

# Maximum theta for which muons can pass through both areas
max_theta = np.arctan(radius / distance_between_areas)

# Rate function considering only allowed angles
def rate_function(p, theta):
    if theta > max_theta:
        return 0  # Muons outside the allowed angular range don't contribute
    return (
        muon_flux(p, theta) / (p**3)  # Momentum normalization
        * effective_area(theta)  # Angular-dependent area
        * np.sin(theta)  # Solid angle integration factor
    )

# Riemann Integration over momentum and angle
pgrid = np.linspace(0.1, 10000, 1000)  # Momentum range
p_step = pgrid[1] - pgrid[0]
costheta_values = np.linspace(1, np.cos(max_theta), 100)  # Range of cos(theta)
costheta_step = (1 - np.cos(max_theta)) / 100
muon_rate = 0 #initializing the muon rate

for costheta in costheta_values: #Riemann Integration over theta and p
    theta = np.arccos(costheta)
    p_integrand = rate_function(pgrid, theta)
    p_integral = np.sum(p_integrand * p_step)
    muon_rate += p_integral * costheta_step

# Multiply by phi symmetry factor
muon_rate *= 2 * np.pi
print("The muon rate for two overlapping circular areas is:", muon_rate, "1/s")

#Plotting the muon rate as a function of vertical separation
# Function to compute muon rate for a given distance
def compute_muon_rate(distance):
    max_theta = np.arctan(radius / distance)  # Max angle for overlapping muons
    costheta_values = np.linspace(1, np.cos(max_theta), 100)  # Range of cos(theta)
    costheta_step = (1 - np.cos(max_theta)) / 100
    muon_rate = 0 #initializing muon rate

    for costheta in costheta_values:
        theta = np.arccos(costheta)
        if theta > max_theta:
            continue  # Skip angles that don't hit both areas
        p_integrand = muon_flux(pgrid, theta) / (pgrid**3) * effective_area(theta) * np.sin(theta)
        p_integral = np.sum(p_integrand * p_step)
        muon_rate += p_integral * costheta_step

    return muon_rate * 2 * np.pi  # Multiply by phi symmetry factor

# Vary the distance between the two areas
distances = np.linspace(5, 100, 20)  # Distances in cm
muon_rates = [compute_muon_rate(d) for d in distances]

# Compute and mark the mDOM reference point
mdom_distance = 41.0  # cm
mdom_rate = compute_muon_rate(mdom_distance)

# Plotting results
plt.figure(figsize=(8, 6))
plt.plot(distances, muon_rates, marker='o', linestyle='-', color='b', label="Muon rate vs. distance")

# Add mDOM reference lines and point
plt.axvline(mdom_distance, color='red', linestyle='--', label=f'mDOM distance = {mdom_distance} cm')
plt.axhline(mdom_rate, color='green', linestyle='--', label=f'mDOM rate ≈ {mdom_rate:.4f} 1/s')
plt.plot(mdom_distance, mdom_rate, 'ro', label='mDOM rate point')

plt.xlabel("Distance between areas (cm)")
plt.ylabel("Muon rate (1/s)")
plt.title("Muon rate through two overlapping areas vs vertical displacement")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

"""The same can be done with the D-Egg module, another module of the Upgrade consisting of only two vertically separated Hamamatsu
R5912-100 PMTs
"""

# D-Egg Constants
diameter = 20.2  # cm
radius = 10.1  # cm
area_circle = np.pi * radius ** 2  # cm^2
distance_between_areas = 53.4  # cm

# Initialize daemonflux
daemonflux = Flux(location='generic')

# Muon flux function with theta dependency
def muon_flux(p, theta):
    return daemonflux.flux(p, theta * 180 / np.pi, 'muflux')  # Convert theta to degrees

# Effective area function considering angle theta
def effective_area(theta):
    return area_circle * np.cos(theta)  # Projection of area due to angle

# Maximum theta for which muons can pass through both areas
max_theta = np.arctan(radius / distance_between_areas)

# Rate function considering only allowed angles
def rate_function(p, theta):
    if theta > max_theta:
        return 0  # Muons outside the allowed angular range don't contribute
    return (
        muon_flux(p, theta) / (p**3)  # Momentum normalization
        * effective_area(theta)  # Angular-dependent area
        * np.sin(theta)  # Solid angle integration factor
    )

# Integration over momentum and angle
pgrid = np.linspace(0.1, 10000, 1000)  # Momentum range
p_step = pgrid[1] - pgrid[0]
costheta_values = np.linspace(1, np.cos(max_theta), 100)  # Range of cos(theta)
costheta_step = (1 - np.cos(max_theta)) / 100
muon_rate = 0

for costheta in costheta_values:
    theta = np.arccos(costheta)
    p_integrand = rate_function(pgrid, theta)
    p_integral = np.sum(p_integrand * p_step)
    muon_rate += p_integral * costheta_step

# Multiply by phi symmetry factor
muon_rate *= 2 * np.pi
print("The muon rate for two overlapping circular areas is:", muon_rate, "1/s")

# Function to compute muon rate for a given distance
def compute_muon_rate(distance):
    max_theta = np.arctan(radius / distance)  # Max angle for overlapping muons
    costheta_values = np.linspace(1, np.cos(max_theta), 100)  # Range of cos(theta)
    costheta_step = (1 - np.cos(max_theta)) / 100
    muon_rate = 0

    for costheta in costheta_values:
        theta = np.arccos(costheta)
        if theta > max_theta:
            continue  # Skip angles beyond the allowed range
        p_integrand = muon_flux(pgrid, theta) / (pgrid**3) * effective_area(theta) * np.sin(theta)
        p_integral = np.sum(p_integrand * p_step)
        muon_rate += p_integral * costheta_step

    return muon_rate * 2 * np.pi  # Multiply by phi symmetry factor

# Vary the distance between the two areas
distances = np.linspace(5, 100, 20)  # Distances in cm
muon_rates = [compute_muon_rate(d) for d in distances]

# Compute and mark the D-Egg reference point
DEgg_distance = 53.4  # cm
DEgg_rate = compute_muon_rate(mdom_distance)

plt.figure(figsize=(8, 6))
plt.plot(distances, muon_rates, marker='o', linestyle='-', color='b', label="Muon rate vs. distance")

# Add D-Egg reference lines and point
plt.axvline(DEgg_distance, color='red', linestyle='--', label=f'D-Egg distance = {mdom_distance} cm')
plt.axhline(DEgg_rate, color='green', linestyle='--', label=f'D-Egg rate ≈ {mdom_rate:.4f} 1/s')
plt.plot(DEgg_distance, mdom_rate, 'ro', label='D-Egg rate point')

plt.xlabel("Distance between areas (cm)")
plt.ylabel("Muon rate (1/s)")
plt.title("Muon rate through two overlapping areas vs vertical displacement")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

"""# Rate through four superposed areas

Using the same method, we can now plot the rate of muons passing through four superposed areas. This is representative of the rate of muons passing through two mDOMs in antarctic ice stacked vertically on top of each other in one string, with their four PMTs aligned
"""

from scipy.interpolate import interp1d #We use this interpolation to calculate the muon rate at the chosen vertical separation that corresponds to mDOM separation in ice

# Constants
diameter = 8.0  # cm
radius = diameter / 2
area_circle = np.pi * radius ** 2  # cm^2
pmt_spacing = 41.0  # cm between PMTs in same mDOM

# Initialize daemonflux
flux_model = Flux(location='generic')

# Muon flux function
def muon_flux(p, theta):
    return flux_model.flux(p, theta * 180 / np.pi, 'muflux')  # degrees

# Effective area projection
def effective_area(theta):
    return area_circle * np.cos(theta)

# Muon Rate function
def rate_function(p, theta, max_theta):
    if theta > max_theta:
        return 0
    return (
        muon_flux(p, theta) / (p**3)
        * effective_area(theta)
        * np.sin(theta)
    )

# distance between mDOMs grid (between two middle PMTs of the stack)
mdom_distances = np.linspace(50, 600, 30)  # cm
muon_rates = []

# Momentum grid
pgrid = np.linspace(0.1, 10000, 1000) #GeV
p_step = pgrid[1] - pgrid[0]

# Riemann Integration over momentum and angle to obtain muon rate for each distance
for d_mdom in mdom_distances:
    total_distance = 2 * pmt_spacing + d_mdom  # Total vertical distance in cm
    max_theta = np.arctan(radius / total_distance)

    costheta_values = np.linspace(1, np.cos(max_theta), 100)
    costheta_step = (1 - np.cos(max_theta)) / 100
    muon_rate = 0

    for costheta in costheta_values:
        theta = np.arccos(costheta)
        p_integrand = rate_function(pgrid, theta, max_theta)
        p_integral = np.sum(p_integrand * p_step)
        muon_rate += p_integral * costheta_step

    muon_rate *= 2 * np.pi  # Azimuthal symmetry
    muon_rates.append(muon_rate)

# Interpolate to get the rate at exactly 300 cm
interp_func = interp1d(mdom_distances, muon_rates, kind='linear')
rate_at_300 = float(interp_func(300))

# Plotting the muon rate and displaying the rate at mDOM separation
plt.figure(figsize=(8, 5))
plt.plot(mdom_distances, muon_rates, marker='o', label='Muon Rate')
plt.axvline(300, color='red', linestyle='--', label='mDOM separation = 300 cm')
plt.axhline(rate_at_300, color='green', linestyle='--', label=f'Rate at 300 cm = {rate_at_300:.2e} 1/s')

plt.xlabel("Distance between mDOMs (cm)")
plt.ylabel("Muon rate through 4 PMTs (1/s)")
plt.title("Muon Rate vs Vertical Distance Between mDOMs")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

