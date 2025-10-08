# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import quad
!pip install daemonflux
import daemonflux
from daemonflux import Flux

"""# Muon rate between two aligned areas

This section builds upon the rate calculated for displaced circles and assumes that they are tilted in a way to be aligned. Once again, the area of the two disks corresponds to the projected circular area of an mDOM. The two modules are located on adjacent strings of the IceCube Upgrade and are vertically displaced one below the other.
"""

#Constants
r = 17.8                     # cm, mDOM radius
area = np.pi * r**2          # cm²
p_min, p_max = 10.0, 1e4     # GeV
n_p, n_theta = 1000, 800     # step size for momentum and angle grid
# Values of mDOM spacing in ice for IceCube Upgrade
x0, z0 = 4500.0, 300.0       # cm

#Momentum Grid
pgrid = np.linspace(p_min, p_max, n_p)
dp    = pgrid[1] - pgrid[0]

#Initialize DAEMONFLUX
flux = Flux(location='generic')

#We first normalize the muon flux and integrate over the momentum grid
def S_of_theta(theta):
    Φ = flux.flux(pgrid, np.degrees(theta), 'muflux') #flux function
    return np.sum(Φ / (pgrid**3)) * dp

#We then define the function that rotates the two disks to be aligned from their initial flat position
def rot_from_ez_to(u): #u is the angle of the rotated circles
    a = np.array([0.0,0.0,1.0])
    b = u/np.linalg.norm(u)
    c = np.dot(a,b)
    if np.isclose(c,1):   return np.eye(3)
    if np.isclose(c,-1):  return -np.eye(3)
    w = np.cross(a,b); s=np.linalg.norm(w)
    W = np.array([[0,-w[2],w[1]],[w[2],0,-w[0]],[-w[1],w[0],0]])
    return np.eye(3)+W+W.dot(W)*((1-c)/(s**2))

#We now can define the rate as a function of horizontal and vertical displacement (x,z),
#accounting for the alignement of the two disks
def rate_two_disks_rotated(x, z):
    #Recomputing disk centers, axis, and rotation depending on (x,z)
    C1 = np.array([0,0,0])
    C2 = np.array([x,0,-z])
    axis = C2 - C1
    n_hat = axis/np.linalg.norm(axis)
    R = rot_from_ez_to(n_hat)

    #Defining the maximum angle that hits both areas
    H = np.linalg.norm(axis)
    theta_max = np.arctan(r/H)

    #Defining the theta grid
    thetas = np.linspace(0, theta_max, n_theta)
    dθ = thetas[1] - thetas[0]

    total = 0.0 #initializing the angular integral
    for θ in thetas:
        #local incoming muons
        dir_local = np.array([np.sin(θ), 0, np.cos(θ)])
        direction = R.dot(dir_local)
        #allowing only muons from disk 1 to disk 2
        if np.dot(direction, n_hat) <= 0:
            continue
        #the intersection condition is named t
        t = np.dot(C2-C1, n_hat) / np.dot(direction, n_hat)
        if t <= 0:
            continue
        hit = C1 + t*direction
        perp = hit - C2
        perp -= np.dot(perp, n_hat)*n_hat
        if np.dot(perp, perp) > r**2:
            continue

        #Integration over all angles that hit both areas
        Sθ = S_of_theta(θ)
        total += Sθ * (area*np.cos(θ)) * (2*np.pi*np.sin(θ)*dθ)

    return total

# compute reference rate
ref_rate = rate_two_disks_rotated(x0, z0)
print(f"Reference rate (x={x0} cm, z={z0} cm): {ref_rate:.2e} 1/s")

#Plot of rate as a function of vertical displacement z

z_vals = np.linspace(50, 1000, 40) #values of vertical displacement
rates_z = [rate_two_disks_rotated(x0, zz) for zz in z_vals] #associated muon rate for each z value with constant x

plt.figure(figsize=(8,5))
plt.plot(z_vals, rates_z, 'b-o', label='Rate vs z')
plt.axvline(z0, color='r', ls='--', label=f'z={z0:.0f} cm')
plt.axhline(ref_rate, color='g', ls='--', label=f'{ref_rate:.2e} 1/s')
plt.xlabel("Vertical separation z (cm)")
plt.ylabel("Muon rate (1/s)")
plt.title("Aligned disks muon rate vs vertical displacement (x = 4500 cm)")
plt.grid(True); plt.legend(); plt.tight_layout(); plt.show()

#Plot of rate as a function of horizontal displacement
x_vals = np.linspace(50, 5000, 40) #values of horizontal displacement
rates_x = [rate_two_disks_rotated(xx, z0) for xx in x_vals] #associated muon rate for each x value with constant z

plt.figure(figsize=(8,5))
plt.plot(x_vals, rates_x, 'm-o', label='Rate vs x')
plt.axvline(x0, color='r', ls='--', label=f'x={x0:.0f} cm')
plt.axhline(ref_rate, color='g', ls='--', label=f'{ref_rate:.2e} 1/s')
plt.xlabel("Horizontal separation x (cm)")
plt.ylabel("Muon rate (1/s)")
plt.title("Aligned disks muon rate vs horizontal displacement (z = 300 cm)")
plt.grid(True); plt.legend(); plt.tight_layout(); plt.show()

"""**Note:** These plots might take a long time to compile (several minutes), as each point individually calculates the rate and realigns the two disks for the given separation.

# Muon Rate through four aligned areas

In analogy to the previous section, we can now calculate the muon rate through four aligned disks. This would represent the rare scenario of a muon traversing the aligned PMTs of two separate mDOMs situated on adjacent strings. This allows to calculate the rate as a function of distance between the two middle PMTs to be compared to how it would behave if the mDOMs were on the same string.
"""

#Constants
r = 4.0                   # cm, PMT radius
A = np.pi * r**2           # cm²
pmt_sep = 41.0             # cm

# IceCube Upgrade spacings
dz_m = 3.0                 # m (vertical separation between mDOMs)
dx_m = 45.0                # m (horizontal string separation)
dz   = dz_m * 100          # cm
dx   = dx_m * 100          # cm

#Momentum grid
pgrid = np.linspace(10.0, 1e4, 1000)
dp    = pgrid[1] - pgrid[0]
#Number of steps for angular grid
n_theta = 500

#Initialize DAEMONFLUX
flux_model = Flux(location='generic')

#We first normalize the muon flux and integrate over the momentum grid
def S_of_theta(theta):
    Φ = flux_model.flux(pgrid, np.degrees(theta), 'muflux') #flux function
    return np.sum(Φ / (pgrid**3)) * dp

#We then define the function that rotates the four disks two by two to be aligned from their initial flat position
def rotation_ez_to(u): #u is the angle of the rotated circles
    a = np.array([0.0,0.0,1.0])
    b = u/np.linalg.norm(u)
    c = np.dot(a,b)
    if np.isclose(c, 1):   return np.eye(3)
    if np.isclose(c,-1):   return -np.eye(3)
    w = np.cross(a,b); s = np.linalg.norm(w)
    W = np.array([[   0,-w[2], w[1]],
                  [ w[2],   0,-w[0]],
                  [-w[1], w[0],   0]])
    return np.eye(3) + W + W.dot(W)*((1-c)/(s**2))

#We now can define the rate as a function of horizontal and vertical displacement (x,z),
#accounting for the alignement of all four disks
def rate_4PMT_stack(dx_cm, dz_cm):
    #Recomputing disk centers, axis, and rotation depending on (x,z)
    C1 = np.array([0.0,       0.0,  +pmt_sep/2])
    C4 = np.array([dx_cm,     0.0, -(dz_cm + pmt_sep/2)])
    v_hat = C4 - C1
    v_hat /= np.linalg.norm(v_hat)
    R = rotation_ez_to(v_hat)

    #Defining the maximum angle that hits both areas
    H = np.linalg.norm(C4 - C1)
    theta_max = np.arctan(r/H)

    #Defining the theta grid
    thetas = np.linspace(0, theta_max, n_theta)
    dθ = thetas[1] - thetas[0]

    total = 0.0 #initializing the angular integral
    for θ in thetas:
        #local incoming muons
        dir_local = np.array([np.sin(θ), 0, np.cos(θ)])
        direction = R.dot(dir_local)

        #allowing only muons that pass from disk 1 to disk 4
        if np.dot(direction, v_hat) <= 0:
            continue

        #the intersection condition is named t
        t = np.dot(C4 - C1, v_hat) / np.dot(direction, v_hat)
        if t <= 0:
            continue
        hit = C1 + t*direction
        perp = hit - C4
        perp -= np.dot(perp, v_hat) * v_hat
        if np.dot(perp, perp) > r**2:
            continue

        #Integration over all angles that hit all areas
        Sθ    = S_of_theta(θ)
        Aproj = A * np.cos(θ)
        dΩ    = 2*np.pi * np.sin(θ) * dθ
        total += Sθ * Aproj * dΩ

    return total

rate_ref = rate_4PMT_stack(dx, dz)
print(f"4-PMT stack reference rate: {rate_ref:.2e} 1/s (dx={dx_m} m, dz={dz_m} m)")

#Plot of muon rate as a function of horizontal displacement x
x_vals_m = np.linspace(0, 50, 20)               # values of horizontal displacement x
rates_x  = [rate_4PMT_stack(xm*100, dz) for xm in x_vals_m] #associated rates for each x value

plt.figure(figsize=(8,5))
plt.plot(x_vals_m, rates_x, 'm-o', label='Rate vs string sep')
plt.axvline(dx_m, color='r', linestyle='--', label=f'x = {dx_m:.0f} m')
plt.axhline(rate_ref, color='g', linestyle='--', label=f'{rate_ref:.2e} 1/s')
plt.xlabel("String separation (m)")
plt.ylabel("Muon rate (1/s)")
plt.title("4-PMT Stack Rate vs horizontal displacement (z = 300 cm)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

#Plot of muon rate as a function of diagonal separation between middle PMTs
H0_m = np.hypot(dx, dz) / 100.0  # reference in meters
H_vals_m = np.linspace(dx_m, H0_m + 5.0, 20)  # values of diagonal displacement
rates_H = []
for Hm in H_vals_m:
    dz_cm = np.sqrt(Hm**2*10000 - dx**2)
    rates_H.append(rate_4PMT_stack(dx, dz_cm)) #associated rates for each H value

plt.figure(figsize=(8,5))
plt.plot(H_vals_m, rates_H, 'b-o', label='Rate vs mid‐PMT H')
plt.axvline(H0_m, color='r', linestyle='--', label=f'H₀ = {H0_m:.2f} m')
plt.axhline(rate_ref, color='g', linestyle='--', label=f'{rate_ref:.2e} 1/s')
plt.xlabel("Middle‐PMT separation (slant H) (m)")
plt.ylabel("Muon rate (1/s)")
plt.title("4-PMT Stack Rate vs diagonal displacement")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

