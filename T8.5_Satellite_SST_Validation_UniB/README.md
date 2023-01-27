[comment]: #![](https://immerse-ocean.eu/img/headers/immerse-header-logo.png) 


<div class="row" style="text-align: center;">
  <div class="column">
   <img src="https://immerse-ocean.eu/img/headers/immerse-header-logo.png"
     alt="Immerse logo"
     style="width:50%" />
  </div>
  <div class="column">
    <img src="https://www.unibo.it/it/immagini/1_UNIBO_Ateneo_vert_pos.jpg/@@images/51c830e4-97ca-4516-bdf8-64ea327bdab3.jpeg"
     alt="PML logo"
     style="width:30%" />
   </div>
</div>

<center><h1>Validation of the NEMO model simulations using Satellite Sea Surface Temperature (SST) measurements</h1></center>


## Introduction
<div style="text-align: justify">
Oil spills in marine environments can have widespread impact and long-term consequences on wildlife, fisheries and coastal habitats. Here, we investigate the impact of fine (submesoscale) scale oceanic motion on the drift of oil. We test whether explicitly representing submesoscale flows allows the more accurate prediction of the pollutant’s advection-diffusion and, therefore, its impacts on the coastal environment. 
</div>
 
## Developments during IMMERSE
<div style="text-align: justify">
As part of the <a href="https://immerse-ocean.eu/">IMMERSE</a> project and <a href="https://marine.copernicus.eu/">CMEMS</a> evolution, a double-nested model experiment have been developed at <a href="https://www.unibo.it/">UNIBO</a> in the offshore waters of the Spanish northwest coast to show the importance of downscaling in the simulation of oil drift.
The <a href="https://www.surf-platform.org/">SURF</a> platform (Trotta et al. 2016), based on the finite differences hydrodynamic <a href="https://www.nemo-ocean.eu/">NEMO</a> code (Madec, 2016; NEMO v3.6 and v4.2) is used to downscale from the parent, comparatively coarse resolution, model to a very high-resolution submesoscale permitting model.
For each nest, the grid spacing decreases by a factor of 4. We use CMEMS reanalysis fields
(<a href="https://data.marine.copernicus.eu/product/IBI_MULTIYEAR_PHY_005_002/description">IBI_REANALYSIS_PHYS_005_002</a>
product) with 1/12° horizontal resolution as the parent model. The NEST1 model domain covers an area of approximately 552 km in longitude by 447 km in latitude, extending from 13.0° W to 6.0° W and from 41.0° N to 45.0° N, with a resolution of 1/48° (∼1638x2316 m). The NEST2 domain extends approximately 200 x 177 km from 11.0° W to 8.5° W and from 42.312° N to 43.9° N, with a resolution of 1/192° (∼417x579 m).
We also consider the impacts of upgrading the NEMO code base to that developed in IMMERSE, v4.2.
</div>

## Purpose of tutorial
<div style="text-align: justify">
The objective of this exercise is to use satellite Sea Surface Temperature (SST) data (<a href="https://data.marine.copernicus.eu/product/SST_MED_SST_L4_NRT_OBSERVATIONS_010_004">SST_MED_SST_L4_NRT_OBSERVATIONS_010_004</a> product) to perform the validation of the high-resolution nested model implemented in the Galicia coastal area (NW of Spain). The simulation started on 7 November 2002 at 00:00 and ran until 17 November 2002 at 24:00. This tutorial consists of one jupyter notebooks (to be found in Downstream-Users-Toolbox/T8.5_Satellite_SST_Validation_UniB/satsst_valid.ipynb) which will:
<ul>
  <li>download and extract SST data from zenodo repository</li>
  <li>read model and satellite SST data</li>
  <li>interpolate the model SST on the satellite SST grid </li>
  <li>generate model and satellite SST maps</li>
  <li>compute and plot the RMSE between measured and interpolated model data</li>
</ul>
</div>


## Install software
<div style="text-align: justify">
<ol>
  <li>
      Install miniconda if needed following https://conda.io/docs/user-guide/install/
Miniconda Distribution includes conda as package manager, Python, the packages they depend on, and a small number of other useful packages.
  </li>
  <li>
      Create conda environment with additional packages needed for tutorial. The following Python packages are needed for running the exercises:
      <div class="row" style="text-align: left;">
      
   | Module name | Description |
| :---: | :---|
| **os** | [Miscellaneous operating system interfaces](https://docs.python.org/3.7/library/os.html) for managing paths, 
| **glob** | [Unix style pathname pattern expansion](https://docs.python.org/3.7/library/glob.html) this module finds all the pathnames matching a specified pattern according to the rules used by the Unix shell |
| **math** | [math](https://docs.python.org/3/library/math.html) is a built-in module in the Python 3 standard library that provides standard mathematical constants and functions.  |
| **numpy** | [NumPy](https://numpy.org/) is the fundamental package for scientific computing with Python and for managing ND-arrays |
| **xarray** | [Xarray](http://xarray.pydata.org/en/stable/) introduces labels in the form of dimensions, coordinates and attributes on top of raw NumPy-like arrays, which allows for a more intuitive, more concise, and less error-prone developer experience. |
| **scipy** |[SciPy](https://docs.scipy.org/doc/scipy/index.html) is an open-source software for mathematics, science, and engineering. It includes modules for statistics, optimization, integration, linear algebra, and more. |
| **matplotlib** |[Matplotlib](https://matplotlib.org/) is a Python 2D plotting library which produces publication quality figures |
| **cartopy** |[Cartopy](https://scitools.org.uk/cartopy/docs/latest/) is a Python package designed for geospatial data processing in order to produce maps and other geospatial data analyses. |
| **python-wget** |[python-wget] is a Python package designed to download remote files from the internet.
   
   You can check all the python packages you have by doing:
<code>
  conda list
</code>
   If some packages are missing you can install them into the current environment using the "conda install" command type the following:
<code>
  conda install numpy xarray matplotlib cartopy python-wget scipy
</code>
   </div>
  </li>
  <li>
   Clone the repository in your local system and navigate to UniB folder
<code>
  git clone https://github.com/immerse-project/Downstream-Users-Toolbox.git
  cd Downstream-Users-Toolbox/T8.5_Satellite_SST_Validation_UniB/
</code>
    </li>
  <li>
  Execute the tutorial
<code>
  jupyter notebook satsst_valid.ipynb
</code>
    </li>
</ol>
   
</div>
