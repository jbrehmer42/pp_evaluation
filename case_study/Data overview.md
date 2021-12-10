# Overview #

## Model output files (forecasts)

These xz-compressed files contain the output of the forecasting
models (see Section 5.2 in the manuscript) in the time period 
from April 16, 2005 to May 26, 2020. The  model outputs are 
matrices. The rows represent model runs (see time stamp file
below) and the columns represent grid cells (see grid cell file
below).

## Earthquake catalog (observations)

Earthquakes in the vicinity of Italy during the testing period.
Each row gives the time, location, depth, and magnitude of one
earthquake.

## Time stamps file

Contain the time stamps corresponding to the model output files,
i.e. entry i in this file specifies the time when the model
output (see model output files above) in row i was produced.
There are a few days where multiple model outputs are available.
On such days only the first model output is used.

## Grid cell file

Locations of grid cells corresponding to the model outputs files,
i.e. entry j in this file specifies the grid cell for which the 
model output (see model output files above) in column j was
issued. Each entry gives the longitude and latitude of the center
of a grid cell and a consecutive number. The edge length of all 
grid cells is 0.1 degrees, so the centers full specify the cells.


