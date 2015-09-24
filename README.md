# art-predictive-model
Predictive Model for Clinical Target Volume using the Kalman Filter

Given an input of a time series T of dth dimensional points P, the PredictiveModel function performs linear quadratic estimation (Kalman filtering) on each point, one dimension at a time, through all samples. Assumes input arguments maintains row correspondence for each point. MATLAB implementation automatically handles arbitrary input sizes wheras the Python implementation needs some minor code changes to modify input sizes.
Outputs: Set of points from final Kalman filter recursion and for MATLAB implementation, raw_data cell array containing all filter outputs

Mathematically:

![Mathematics](https://raw.githubusercontent.com/yazanobeidi/art-predictive-model/master/formulas.png?raw=true)

MATLAB:
Runtime example can be executed with runPredictiveModel script from command.

PYTHON:
Place plan.roi and predModel.py in the same directory. Run predModel.py. See output.md or newartplan.roi for output. Also see comments in predModel.py to output to console. 

This algorithm was designed as part of adaptive radiation therapy for clinical
target volume determination. 

Copyright 2015 Yazan Obeidi
