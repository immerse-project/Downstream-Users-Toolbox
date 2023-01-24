[comment]: #![](https://immerse-ocean.eu/img/headers/immerse-header-logo.png) 

<div class="row">
  <div class="column">
   <img src="https://immerse-ocean.eu/img/headers/immerse-header-logo.png"
     alt="Immerse logo"
     style="width:30%" />
  </div>
  <div class="column">
    <img src="https://www.hereon.de/cms60/res/assets/logos/hereon_logo.svg"
     alt="hereon logo"
     style="width:30%" />
   </div>
</div>

<center><h1> Validation of structured grid (NEMO) and unstructured grid (SCHISM) model simulations using lagrangian FerryBox measurements
  
  
## Introduction
The goal of this tutorial is to validate surface temperature and salinity outputs from different models using lagrangian data from a FerryBox, travelling from the outer estuaty towards the open ocean. Within the IMMERSE project, the approach has been used to evaluate the capabability of a newly developed high resolution (400m) NEMO configuration for the Southern North Sea to replicate the outer estaury dynamics,
compared to the coarser CMEMS operational prdouct AMM15, and the locally higher resolving pre-operational SCHISM based setup from Hereon.
  
### Background
Ocean models are valuable tools to predict future (forecast) or reconstruct past (hindcast) events, gain undestranding of physical processes, and test 'what if scenairos' in the context of adaptation and management strategies. Credibility regarding the reasonable functioning of such models is gained from the validation against observational data, indicating how good the model can replicate the observed ocean state.
  
  
### Developments during IMMERSE
In the context of IMMERSE at Hereon a NEMO configuration for Southern North Sea (SNS) with a high resolution of 400m was setup, based on the new NEMO.v4.2_RC code introducing the feature of wetting and drying. 
This setup was tested and validated in comparison to other models, and later on coupled with the WAM wave model.
  
### Purpose of tutorial
This tutorial shows how to plot maps of different variables from the different models with an overlay of the measurement data, perform interpolations to colocate the model outputs with the FerryBox in time an space, and compute and plot some basic error statistics.

The repository includes one jupyter notebook (ferrybox_validation.ipynb) and a python script within in a subfolder lib, containing functions to facilitate working with the unstructured grid SCHISM output (schism.py)).

## Install software
Clone the repository in your local system and navigate to ferrybox_validation folder

```
git clone https://github.com/immerse-project/Downstream-Users-Toolbox.git
cd Downstream-Users-Toolbox/T8.2_Ferrybox_Validation_Hereon
```

## Get the Data

```
cd Downstream-Users-Toolbox/T8.2_Ferrybox_Validation_Hereon
wget https://zenodo.org/record/7521488/files/data.zip
unzip data.zip
```
## Execute the tutorial

```
jupyter notebook ferrybox_validation.ipynb
```



