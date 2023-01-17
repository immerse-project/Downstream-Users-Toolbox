[comment]: #![](https://immerse-ocean.eu/img/headers/immerse-header-logo.png) 

<div style=\max-width:50%; display:inline-block;\ \>
 </div>
 
<img src="https://immerse-ocean.eu/img/headers/immerse-header-logo.png"
     alt="Immerse logo"
     style="float: left;text-align: center; width: 350px; margin-right: 10px;" />
<img src="https://www.pml.ac.uk/Content/Images/PMLLogo.svg"
     alt="PML logo"
     style="float: right;text-align: right; width: 350px; margin-right: 10px;" />

<center><h1> Evaluation and Analysis of estuarine hydrodynamics with the unstructured grid model FVCOM

## Introduction
The goal of this tutorial is to analyze estuarine dynamics within the Total Exchange Flow (TEF) framework as orginally developed by McCready (2011) and subsequently extended in Wang et al. (2017) and Burchard et al. (2019). The TEF extends the work of Knudsen on estuarine salinity volume exchanges to discretised such exchange in the salinity space. By doing so, one can derive both the advective and diffusive fluxes as well as quantify turbulent mixing through changes in the salinity variance. This framework can be used to validate model performance or to assess changes in model setup. Within the IMMERSE project, the approach has been used to evaluate the effects of using different CMEMS products (AMM15 vs AMM7 NorthWest Shelf operational products) on the dynamics of the Tamar estuary (UK) and their subsequent impact on biogeochemistry using the unstructured grid hydrodynamic model FVCOM (Chen et al 2003, Sims et al., 2022) and the biogeochemical ERSEM model (Butenschön et al., 2016).

### Background
Mixing of freshwater and ocean waters is an important process that occurs in the coastal region and very prominently in estuarine environments. With estuaries concentrating a large percentage of the world population, they recieve large volumes of pollutants resulting from anthropogenic activities associated with land use such as organic and inorganic particulates (from soil erosion caused by deforestation, overgrazing, and other poor farming practices) or excessive nutrients (from fertilizer and animal  wastes including aquaculture) and with high population density such as nutrients (seawage), pollutants including heavy metals (industrial waste and urban runoff) or chlorinated hydrocarbons (e.g. from insecticides). 

The description of how these pollutants mix within estuaries is therefore crucial to understanding the potential impacts on coastal and nearshore ecosystems. The use of bespoke high resolution models is one approach to studying this problem. Because they require substantial computing resources, they cover small coastal sections and are forced at their ocean boundaries by larger coarser models. 

### Developments during IMMERSE
As part of the IMMERSE project and CMEMS evolution, the North West European shelf operational model increased its horizontal resolution from 7km (AMM7) to 1.5 km (AMM15) and additional variables were made available at hourly frequency instead of daily. FVCOM model code was enhanced to perform biogeochemical simulations using previously stored hydrodynamics results. 

### Purpose of tutorial
This tutorial shows how to define transects and extract the relevant variables from FVCOM model outputs. It then interpolates all variables to a common grid and calculates the subtidal, time-varying normal volume flux in each transect for each salinity bin. It then displays the residual fluxes and the mean flux per salinity bin accross each transect. Finally it shows the variability in mixing as depicted by the salinity variance metric. 

The repository includes one jupyter notebook (mixing_calc.ipynb) and two python scripts with functions to facilitate working with FVCOM output (aux_func_nopyfvcom.py) and for filtering the tides from the model output (tidal_filter.py). A short model output file is accessed through thredds to perform the calculations step by step. Pre-calculated quantities on a longer timeseries is included here (cross_sect_defn.npy).

## Install software
Clone the repository in your local system and navigate to PML folder

```
git clone https://github.com/immerse-project/Downstream-Users-Toolbox.git
cd Downstream-Users-Toolbox/PML
```

## Execute the tutorial

```
jupyter notebook mixing_calc.ipynb
```

## References 
Butenschön, M., Clark, J., Aldridge, J.N., Allen, J.I., Artioli, Y., Blackford, J., Bruggeman, J., Cazenave, P., Ciavatta, S., Kay, S., Lessin, G., van Leeuwen, S., van der Molen, J., de Mora, L., Polimene, L., Sailley, S., Stephens, N., Torres, R., 2016. ERSEM 15.06: a generic model for marine biogeochemistry and the ecosystem dynamics of the lower trophic levels. Geoscientific Model Development 9, 1293–1339. https://doi.org/10.5194/gmd-9-1293-2016

Burchard, H., Lange, X., Klingbeil, K., & MacCready, P. 2019. Mixing Estimates for Estuaries, Journal of Physical Oceanography, 49(2), 631-648. Retrieved Jan 12, 2023, from https://journals.ametsoc.org/view/journals/phoc/49/2/jpo-d-18-0147.1.xml

Chen, C., Liu, H., Beardsley, R.C., 2003. An Unstructured Grid, Finite-Volume, Three-Dimensional, Primitive Equations Ocean Model: Application to Coastal Ocean and Estuaries. Journal of Atmospheric and Oceanic Technology 20, 159–186. https://doi.org/10.1175/1520-0426(2003)020%3C0159:AUGFVT%3E2.0.CO;2

MacCready, P., 2011: Calculating estuarine exchange flow using isohaline coordinates. J. Phys. Oceanogr., 41, 1116–1124, https://doi.org/10.1175/2011JPO4517.1.

Sims, R.P., Bedington, M., Schuster, U., Watson, A.J., Kitidis, V., Torres, R., Findlay, H.S., Fishwick, J.R., Brown, I., Bell, T.G., 2022. Tidal mixing of estuarine and coastal waters in the western English Channel is a control on spatial and temporal variability in seawater CO2. Biogeosciences 19, 1657–1674. https://doi.org/10.5194/bg-19-1657-2022


