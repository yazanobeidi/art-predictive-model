function [ model, raw_data ] = predictiveModel( varargin )
%PREDICTIVEMODEL Perform Kalman Filter on a time series of nth dimensional points
%   Given a variable input of a time series set of nth dimensional points, this
%   function performs a Kalman Filter on each point, one dimension at a time.
% 	Currently, written to handle up to 11 input samples, more can be added with ease.
%   Assumes input argument maintains row correspondence for each point.
%   Outputs: model n x 3 matrix containing last recursive Kalman Filter
%   iteration and raw_data n x 3 cell array containing all Filter outputs
%
%   Logic:
%   let Pn = (xn,yn,zn) for all n > 0
%   let Fm = [Pn] for all m > 0
%   Fm = F(m-1) +- delta 
%   Define Tn(Pn,Dd) = [F0(P0,D0);... Fm(Pn,Dn)] for all DdCPn, d = (1,2,3)
%   Algorithm:
%   output_n_d = kalman(Tn(Pn,Dd))
%   model(n,d) = output_n_d(end) where model is a set of n points, for all n, d > 0

% Get input arguments, can continue with pattern to add additional input arguments
if nargin >= 1
    f0 = varargin{1};
if nargin >= 2
    f1 = varargin{2};
if nargin >= 3
    f2 = varargin{3};
if nargin >= 4
    f3 = varargin{4};
if nargin >= 5
    f4 = varargin{5};
if nargin >= 6
    w1 = varargin{6};
if nargin >= 7
    w2 = varargin{7};
if nargin >= 8
    w3 = varargin{8};
if nargin >= 9
    w4 = varargin{9};
if nargin >= 10
    w5 = varargin{10};
if nargin >= 11
    w6 = varargin{11};
end
end
end
end
end
end
end
end
end
end
end

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
        switch nargin
            case 1
                Tn_a = f0(Pn,Dd);
            case 2
                Tn_a = [f0(Pn,Dd);f1(Pn,Dd)];
            case 3
                Tn_a = [f0(Pn,Dd);f1(Pn,Dd);f2(Pn,Dd)];
            case 4
                Tn_a = [f0(Pn,Dd);f1(Pn,Dd);f2(Pn,Dd);f3(Pn,Dd)];
            case 5
                Tn_a = [f0(Pn,Dd);f1(Pn,Dd);f2(Pn,Dd);f3(Pn,Dd);f4(Pn,Dd)];
            case 6
                Tn_a = [f0(Pn,Dd);f1(Pn,Dd);f2(Pn,Dd);f3(Pn,Dd); ...
                    f4(Pn,Dd);w1(Pn,Dd)];
            case 7
                Tn_a = [f0(Pn,Dd);f1(Pn,Dd);f2(Pn,Dd);f3(Pn,Dd); ...
                    f4(Pn,Dd);w1(Pn,Dd);w2(Pn,Dd)];
            case 8
                Tn_a = [f0(Pn,Dd);f1(Pn,Dd);f2(Pn,Dd);f3(Pn,Dd); ...
                    f4(Pn,Dd);w1(Pn,Dd);w2(Pn,Dd);w3(Pn,Dd)];
            case 9
                Tn_a = [f0(Pn,Dd);f1(Pn,Dd);f2(Pn,Dd);f3(Pn,Dd); ...
                    f4(Pn,Dd);w1(Pn,Dd);w2(Pn,Dd);w3(Pn,Dd); ...
                    w4(Pn,Dd)];
            case 10
                Tn_a = [f0(Pn,Dd);f1(Pn,Dd);f2(Pn,Dd);f3(Pn,Dd); ...
                    f4(Pn,Dd);w1(Pn,Dd);w2(Pn,Dd);w3(Pn,Dd); ...
                    w4(Pn,Dd);w5(Pn,Dd)];
            case 11
                Tn_a = [f0(Pn,Dd);f1(Pn,Dd);f2(Pn,Dd);f3(Pn,Dd); ...
                    f4(Pn,Dd);w1(Pn,Dd);w2(Pn,Dd);w3(Pn,Dd); ...
                    w4(Pn,Dd);w5(Pn,Dd);w6(Pn,Dd)];
        end % can continue with pattern to add additional input arguments
        
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

% Written by Yazan Obeidi July 2015, Princess Margaret Cancer Centre
