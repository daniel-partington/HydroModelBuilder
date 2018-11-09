[![DOI](https://zenodo.org/badge/63133420.svg)](https://zenodo.org/badge/latestdoi/63133420)

# HydroModelBuilder

HydroModelBuilder is made up of a set of utilities to build, run, and optimise numerical groundwater models.
The primary use case of HydroModelBuilder is to generate groundwater models within an integrated environmental
modelling process.

It is currently under development and supports the generation of MODFLOW models. Future capabilities will expand the
supported groundwater modelling platforms, include surface water models, and the ability generate model emulators.
Collaboration is encouraged and appreciated.

HydroModelBuilder currently consists of two main modules: *GWModelBuilder* and *GWModelManager*

GWModelBuilder
-------------------------
Module to:
* Build a database of model relevant data
  (e.g. rasters of aquifers layers, shapefile of river, time series of river level etc)
* Build a model mesh
* Translate model relevant data to a model mesh
  (e.g. mapping GIS data to model mesh)
* Package model data

GWModelManager
-------------------------
Module to:
* Load and run models
* Optimise models and conduct uncertainty analysis
* In future: Emulate numerical models
