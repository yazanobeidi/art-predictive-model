# art-predictive-model
Predictive Model for Clinical Target Volume using the Kalman Filter in MATLAB

Given a variable input of a time series set of nth dimensional points, this
function performs a Kalman Filter on each point, one dimension at a time.
Currently, written to handle up to 11 input samples, more can be added with ease.
Assumes input argument maintains row correspondence for each point.
Outputs: model n x 3 matrix containing last recursive Kalman Filter
iteration and raw_data n x 3 cell array containing all Filter outputs

Runtime example can be executed with runPredictiveModel script from command.
Example scenario is currently a cube shifting over 5 'samples', will be updated
with something more interesting in the near future.

Logic:
* let Pn = (xn,yn,zn) for all n > 0
* let Fm = [Pn] for all m > 0
* Fm = F(m-1) +- delta
* Define Tn(Pn,Dd) = [F0(P0,D0); F1(P1,D1) ... Fm(Pn,Dn)] for all DdCPn, d > 0
  
Algorithm: {

  output_n_d = kalman(Tn(Pn,Dd)) 

  model(n,d) = output_n_d(end) where model is a set of n points, for all n, d > 0

}
