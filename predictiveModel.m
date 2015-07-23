%	Copyright 2015 Yazan Obeidi
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License along
%   with this program; if not, write to the Free Software Foundation, Inc.,
%   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

function [ model, raw_data ] = predictiveModel( varargin )
%PREDICTIVEMODEL Perform Kalman Filter on a time series of nth dimensional points
%   Given a variable input of a time series set of nth dimensional points, this
%   function performs a Kalman Filter on each point, one dimension at a time.
% 	Can handle any amount of inputs. Assumes input argument maintains row 
%	correspondence for each point. Outputs model n x 3 matrix containing last 
%	recursive Kalman Filter iteration and raw_data n x 3 cell array containing all 
%	Filter outputs.
%
%   Logic:
%   let Pn = (xn,yn,zn) for all n > 0
%   let Fm = [Pn] for all m > 0
%   Fm = F(m-1) +- delta 
%   Define Tn(Pn,Dd) = [F0(P0,D0);... Fm(Pn,Dn)] for all DdCPn, d > 0
%   Algorithm:
%   output_n_d = kalman(Tn(Pn,Dd))
%   model(n,d) = output_n_d(end) where model is a set of n points, for all n, d > 0

% Get number of data points and number of dimensions
n = size(varargin(1)); d = n(2); n = n(1);
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
        Tn_a = zeros(nargin,1);
        for i = 1:nargin % iterate through all inputs
        	temp = varargin(i);
        	var = temp(1,1)) % since matlab does not allow varargin(i,1)
        	Tn_a(i,1) = var(Pn,Dd);
		end
        
        % Case specific Kalman Filter variables
        output_n_d = []; % initialize results store
        xi = Tn_a(1,1); % initial state vector
        
        % Perform Kalman Filter
        output_n_d = kalmanf(Tn_a,A,xi,Pi,B,uk,wk,Q,R,H,output_n_d,t);
        
        % Record results
        model(Pn,Dd) = output_n_d(end);
        raw_data{Pn,Dd} = output_n_d;
    end
end
return;
end

% Written by Yazan Obeidi July 2015 at the Princess Margaret Cancer Centre
