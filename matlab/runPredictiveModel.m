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

%%% Kalman Filter Run Script with Mock Patient Data
% This script applies the kalman filter on a test volume time series
% - as opposed to entering these through command

% Create mock data
% 1x1x1 cube with corner at orgin:
initial_scan = [0,0,0; 1,0,0; 0,1,0; 1,1,0; 0,0,1; 1,0,1; 0,1,1; 1,1,1];
frac_one = initial_scan * .90 + rand(1);
frac_two = frac_one * .80 + rand(1);
frac_three = frac_two * 1.10 - rand(1);
frac_four = frac_three - frac_two + frac_one;

[model, raw_data] = predictiveModel(initial_scan, frac_one, frac_two, frac_three, frac_four);

scatter3(initial_scan(:,1),initial_scan(:,2),initial_scan(:,3),'filled','r')
hold on
scatter3(model(:,1),model(:,2),model(:,3),'filled','b')

% Written by Yazan Obeidi July 2015 at the Princess Margaret Cancer Centre
