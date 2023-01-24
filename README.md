
<div class="row">
  <div class="column">
   <img src="https://immerse-ocean.eu/img/headers/immerse-header-logo.png"
     alt="Immerse logo"
     style="width:30%" />
  </div>

# Downstream User Validation and Analysis Toolbox
Jason Holt, National Oceanography Centre, UK

Joanna Staneva, Sebastian Grayek, Benjamin Jacob, Helmoltz-Zentrum Hereon, Germany

Erik Van Sebille, Siren Rühs, University of Utrecht, Netherlands

Ricardo Torres, Plymouth Marine Laboratory, UK

Francisco Trotta, University of Bologna, Italy

The Copernicus Marine Service (CMS; https://marine.copernicus.eu/) by design focuses on activities that are best delivered at a global to regional scale by a pan-European effort. This naturally leads it to generic products that are widely applicable across a diverse range of sectors. On the other hand, end users generally require bespoke information tailored to the particular question at hand, often crossing discipline boundaries, at finer (national to local) scale and requiring more in-depth analysis or processing. This gap is closed by downstream services that provide an additional layer of modelling and analysis. 


This respository provides a set of model validation and analysis tools in the form of python jupyter Notebooks.
They were developed as part of WP8 in the IMMERSE Horizons 2020 project and offer best practice examples of how to
share analysis and assessment approaches among downstream users. They relate to the four downstream case studies in IMMERSE (see Deliverable D8.2 Outcomes of the Case Studies; https://immerse-ocean.eu/deliverables/):

T8.2 Coastal processes in German Bight (Hereon) 

T8.3 Marine plastic litter transport (University of Utrecht)

T8.4 Water quality modelling of the Tamar Estuary and adjacent coast (Plymouth Marine Laboratory)

T8.5 Pollution transport by Submesoscales in the open ocean (University of Bologna)

Each case study notebook is accompanied by a readme file, providing a complete description of how to setup and run the Notebook,
including installing an appropriate python environment and accessing the data. Data for T8.2, T8.3 and T8.5 are available on zenodo.org, while data for T8.4 in provided by an on-line THREDDS server. 

The case study notebooks can run easily on local computer resources (e.g. windows, linux or mac desktops/laptops), with an suitable python environment installed, e.g. via anaconda (www.anaconda.com) They can also run on cloud based systems, such as the WEkEO Data and Information Access Services (DIAS; https://www.wekeo.eu/) system. Case study Notebooks T8.2, 4, and 5 run on the Machine Learning python environment provided with WEkEO. The T8.3 Notebook includes an application of the oceanparcels Lagrangian model. The required data volume for the accompanying NEMO model data and the required python environment together exceed the storage capacity available for free WEkEO DIAS access.

The Notebooks each provide detailed instructions for running, but the necessary prerequisites to run all of them can 
be installed at once creating a common environment containing all the required python packages via the following:

```
git clone git@github.com:immerse-project/Downstream-Users-Toolbox.git

cd Downstream-Users-Toolbox/
```
Install a common conda environment
```
conda create -n immerse_env -c conda-forge --file conda_install.txt
```
Make this environment visible as kernel for jupyter notebooks
```
conda activate immerse_env
python -m ipykernel install --user --name immerse_env 

```




