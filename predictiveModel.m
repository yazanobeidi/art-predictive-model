function [ model, raw_data ] = predictiveModel( f0, f1, f2, f3, f4 )
%PREDICTIVEMODEL Perform Kalman Filter on 5 time series of 3D points
%   Input initial set of points and 4 other time series evolutions, this
%   function performs a Kalman Filter on each point, dimension at a time.
%   Assumes input argument maintains row correspondence for each point.
%   Outputs model n x 3 matrix containing last recursive Kalman Filter
%   iteration and raw_data n x 3 cell array containing all Filter outputs
%
%   Logic:
%   let Pn = (xn,yn,zn) for all n > 0
%   let Fm = [Pn] for all m > 0
%   Fm = F(m-1) +- delta 
%   Define Tn(Pn,Dd) = [F0(P0,D0);... Fm(Pn,Dn)] for all DdCPn, d = (1,2,3)
%   Algorithm:
%   output_n_d = kalman(Tn(Pn,Dd))
%   model(n,d) = output_n_d(end) where model is a set of n points

% Get number of data points and number of dimensions
n = size(f0); d = n(2); n = n(1);
% Initialize model return variable for speed
model = zeros(n,d);
% Initialize raw data return variable cell array
raw_data = cell(n,d);

% Seed Kalman Filter variables
A = 1; % state transition coefficient matrix
B = 0; % control coefficient matrix
uk = 0; % control signal vector
Pi = 1; % initial covariance matrix
wk = 0; % process noise vector (for each state parameter)
Q = cov(wk); % process covariance
R = 0.1;% measurement covariance
H = 1; % transformation matrix
t = nargin; % set recursive counter.. number of input sets

for Pn = 1:n % For all points
    for Dd = 1:d % For all 3 dimensions

        % Define time series T(Pn,An)
        Tn_a = [f0(Pn,Dd);f1(Pn,Dd);f2(Pn,Dd);f3(Pn,Dd);f4(Pn,Dd)];
        
        % Case specific Kalman Filter variables
        output_n_d = []; % initialize results store
        xi = Tn_a(1,1); % initial state vector
        
        % Perform recursive Kalman Filter
        output_n_d = kalmanf(Tn_a,A,xi,Pi,B,uk,wk,Q,R,H,output_n_d,t);
        
        % Record results
        model(Pn,Dd) = output_n_d(end);
        raw_data{Pn,Dd} = output_n_d;
    end
end
return;
end

% Written by Yazan Obeidi July 2015, Princess Margaret Cancer Centre
