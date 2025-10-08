# Muon-Rate-for-Different-Geometries
This project aims to calculate the rate of muons passing through different simplified geometries meant to simulate the module configurations of the IceCube Upgrade. It serves as a first step for a future full simulation. 

The code in this repository is written in the context of obtaining the bachelor of Science in Physics at the Humboldt University Berlin, Germany, titled "Estimating the Rate of Atmospheric Muons Crossing mDOMs in
the IceCube Upgrade". It uses the resources made available by the Deutsches Elektron-Synchrotron and the IceCube Collaboration. The main package used in all the code to estimate the flux of muons is the  DAta-drivEn
MuOn-calibrated atmospheric Neutrino Flux (DAEMONFLUX), found here: https://github.com/mceq-project/daemonflux . 

This project has multiple iterative steps that each describe different geometries that are useful to simulate for the IceCube Upgrade. However, the way the muon rate is calculated in each remains analogous: 

We first define the constants, i.e. the dimensions of the circular areas and their separation, be it vertical and/or horizontal. 

Then, we initialize the muon flux using DAEMONFLUX and establish a grid size for the angle and the momentum, as the flux is given in (p/(GeV/c))^3/ ((GeV/c) s sr cm^2). 

Using this grid, we Riemann Integrate the normalized flux over the momentum and the solid angle to obtain the muon rate passing through the area. 

This method allows to vary the parameters such as area, distance between areas and displacement for each geometry chosen to accurately model the IceCube Upgrade. 

## Requirements 
* `Python > 3.7`, `numpy`, `scipy`, `daemonflux`, `pandas` 
* `matplotlib` for plots
