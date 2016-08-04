# HydroModelBuilder

Under development.

Utilities to build, run, optimise and emulate numerical surface water-groundwater models

Consists of two main modules: *GWModelBuilder* and *GWModelManager*

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
* Emulate numerical models