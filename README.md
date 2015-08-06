# art-predictive-model
Predictive Model for Clinical Target Volume using the Kalman Filter in MATLAB

Given a variable input of a time series set of nth dimensional points, this
function performs a Kalman Filter on each point, one dimension at a time.
Written to handle any arbitrary number of inputs. Assumes input arguments maintains row correspondence for each point.
Outputs: filtered model n x 3 matrix containing last recursive Kalman Filter
iteration and raw_data n x 3 cell array containing all Filter outputs

Runtime example can be executed with runPredictiveModel script from command.
Example scenario is currently a cube shifting over 5 'samples', will be updated
with something more interesting in the near future. To run, open path in MATLAB and type 'runPredictiveModel'.

This algorithm was designed as part of adaptive radiation therapy for clinical
target volume determination. 

Mathematically:

![Mathematics](https://raw.githubusercontent.com/yazanobeidi/art-predictive-model/master/formulasAndLogic.png?raw=true)

Copyright 2015 Yazan Obeidi
