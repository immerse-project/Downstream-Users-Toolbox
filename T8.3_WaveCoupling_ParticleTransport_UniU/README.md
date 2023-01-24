<div class="row">
  <div class="column">
   <img src="https://immerse-ocean.eu/img/headers/immerse-header-logo.png"
     alt="Immerse logo"
     style="width:30%" />
  </div>
  
# Assessing impact of NEMO wave coupling on surface particle dispersal simulations

## Introduction
The goal of this tutorial is to analyze the impact of waves for surface particle dispersal simulations. Specifically, the impact of new options for the representation of wave-current interactions in NEMO implemented during the IMMERSE project, on particle simulations with the Parcels software [OceanParcels.org; Lange and van Sebille, 2017; Delandmeter and van Sebille 2019] are tested in a case study for the Mediterranean Sea.
### Background
The knowledge of how seawater moves around in the global ocean and transports tracers and particulates, is crucial for solving many outstanding issues in physical oceanography and climate science. Due to limited available observations, seawater pathways are often estimated by evaluating virtual particle trajectories inferred from velocity fields computed with ocean models. The quality of these Lagrangian analyses strongly depends on how well the underlying ocean model represents the ocean circulation features of interest. Under influence of surface waves, particles do not only move according to the Eulerian current velocity, but additionally experience a net drift in the direction of wave propagation, called Stokes drift. Moreover, the presence of waves alters the Eulerian current field, which thus can be divided into a wave-driven component and a non-wave-driven component. Notably, the wave-driven Eulerian currents tend to act in the opposing direction of Stokes drift [e.g., van den Bremer and Breivik 2017; Higgins 2020]. 
However, until now, large-scale Lagrangian dispersal simulations mainly used velocity output from ocean models that do not explicitly resolve surface waves, which implies that the wave impact is either not included or must be approximated. Most commonly, for surface particle dispersion simulations, a superposition of the Eulerian currents from an ocean-only model and Stokes drift from an independently run wave model are used [e.g., Onink et al. 2019; Bosi et al., 2021]. 
### Developments during IMMERSE
As part of the IMMERSE project, the options for the representation of wave-current interactions in the NEMO ocean model have been improved and tested, e.g., in a regional high-resolution (1/24° horizontal resolution) model configuration for the Mediterranean Sea. These tests include a coupled ocean-wave model simulation and a complimentary ocean-only simulation, described in detail in deliverable D5.7 “Assessment of wave-current effects on the circulation in the Med-MFC system” (the simulations make use of the NEMO v4.2-RC ocean model, the Wave Watch 3 v.6.07 wave model, the OASIS3-MCT coupler, and ECMWF atmospheric fields). 
Additionally, the Lagrangian software Parcels (used to calculate particle trajectories from model velocity output) has been adjusted to better interpolate velocities on the ORCA C-grid used in NEMO [Delandmeter and van Sebille 2019]. The new interpolation scheme sustains impermeability at the coastal boundary (while with the previous interpolation schemes designed for A-grids, a non-negligible number of particles would artificially beach due to non-zero velocities from ocean into land grid cells) and has significantly improved the accuracy of Lagrangian simulations in NEMO flow fields on ORCA grids. 
### Purpose of tutorial
This tutorial offers some basic analysis tools to address the following questions: How do waves impact simulated surface particle dispersal, and what is the relative impact of Stokes drift and wave-driven Eulerian currents? How well can the wave impact be approximated by the superposition of Eulerian mean and Stokes drift velocity fields obtained from independently run ocean and wave models?
It consists of two jupyter notebooks:
1. Parcels_CalcTraj.ipynb: Calculates Lagrangian particle trajectories based on velocity output (either only Eulerian currents or Eulerian currents plus Stokes drift) from ocean only as well as coupled ocean-wave model simulation by making use of Parcels
2. CompTraj_uncoupledVScoupled.ipynb: Compares Lagrangian particle trajectories calculated from ocean only and coupled ocean-wave model simulations to assess impact of waves on surface particle dispersal


## Install software
1. Install miniconda if needed following https://conda.io/docs/user-guide/install/
2. Create conda environment with additional packages needed for tutorial (parcels: https://oceanparcels.org/, jupyter: https://jupyter.org/, xhistogram: https://xhistogram.readthedocs.io/en/latest/installation.html)
```
conda activate root
conda create -n py3_parcels-waves -c conda-forge parcels cartopy ffmpeg jupyter xhistogram
```
## Install software
Clone the repository in your local system and navigate to Wave Coupling and Particle trackign folder:

```
git clone https://github.com/immerse-project/Downstream-Users-Toolbox.git
cd Downstream-Users-Toolbox/T8.3_WaveCoupling_ParticleTransport_UniU
```

## Gather example data
```

mkdir data
cd data
mkdir domain
mkdir surface_TKE_UNC 
mkdir surface_TKE_CO
cd domain
### ->  download grid data ###
cd ../surface_TKE_UNC
### -> download surface velocity data from uncoupled simulation (surface_TKE_v42RC) ###
cd ../surface_TKE_CO
### -> download surface velocity data from coupled simulation (surface_TKE_CO_FORCE_MIX_LC015) ###
cd .. #back to data folder
```

## Execute the tutorial
Create folder for trajectory and figure output
```
mkdir LAtrajectories
mkdir ../figures
```
Activate environment and start jupyter notebook
```

conda activate py3_parcels-waves
export CC==gcc # necessary when Parcels is installed on MAC
jupyter notebook
```
Open and execute notebooks, first Parcels_CalcTraj.ipynb, then CompTraj_uncoupledVScoupled.ipynb

##References
Bosi, S., Broström, G., & Roquet, F. (2021). The Role of Stokes Drift in the Dispersal of North Atlantic Surface Marine Debris. Frontiers in Marine Science, 8, 1137

Delandmeter, P., and E. van Sebille (2019), The Parcels v2.0 Lagrangian framework: new field interpolation schemes, Geosci. Model Dev., 12(8), 3571-3584, doi:10.5194/gmd-12-3571-2019.

Higgins, C., Vanneste, J., & van den Bremer, T. S. (2020). Unsteady Ekman-Stokes Dynamics: Implications for Surface Wave-Induced Drift of Floating Marine Litter. Geophysical Research Letters, 47(18). https://doi.org/10.1029/2020GL089189

Lange, M., and E. van Sebille (2017), Parcels v0.9: prototyping a Lagrangian ocean analysis framework for the petascale age, Geosci. Model Dev., 10(11), 4175-4186, doi:10.5194/gmd-10-4175-2017.

Onink, V., Wichmann, D., Delandmeter, P., & van Sebille, E. (2019). The Role of Ekman Currents, Geostrophy, and Stokes Drift in the Accumulation of Floating Microplastic. Journal of Geophysical Research: Oceans, 124(3), 1474–1490. 

van den Bremer, T.S., & Breivik, Ø. (2017). Stokes drift. Phil. Trans. R. Soc. A, 376:20170104. http://dx.doi.org/10.1098/rsta.2017.0104 



