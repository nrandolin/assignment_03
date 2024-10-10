clear all
filterparams.min_xval = 0;
filterparams.max_xval = 1;

%% GLOBAL FORWARD EULER

%% LOCAL FORWARD EULER
t = 0.5;
h_list = linspace(10E-5,10,100);
tspan = [0,t];
X0 = [0; 1];
local_error = [];
difference = [];

for i = 1:length(h_list)
   h_ref = h_list(i);
   % Calculate numperical x value
   [t_list,X_list,h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
   X_numerical = X_list(end);
   % calculate the real x value
   X_true = solution01(t+h_ref);
   local_error = [local_error, abs(X_numerical-X_true)];
   difference = [difference, abs(solution01(t) - X_true)];
end

figure()
% Log-log of local truncation error
loglog(h_list, local_error, '.b', 'MarkerSize', 12)
title("Local Truncation Error, Forward Euler")
xlabel("Step Size")
ylabel("Error")
[p,k] = loglog_fit(h_list,local_error, filterparams);
hold on
% Difference Line
loglog(h_list, difference, '.r', 'MarkerSize', 12)
% Fit Line
loglog(h_list, k*h_list.^p, '--m', 'LineWidth', 1.5)
legend("Local Truncation Error", "Difference", "Fit Line")
%% LOCAL MIDPOINT
t = 0.5;
h_list = linspace(10E-5,10,50);
tspan = [0,t];
X0 = [0; 1];
local_error = [];
difference = [];

for i = 1:length(h_list)
   h_ref = h_list(i);
   % Calculate numperical x value
   [t_list,X_list,h_avg, num_evals] = explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
   X_numerical = X_list(end);
   % calculate the real x value
   X_true = solution01(t+h_ref);
   local_error = [local_error, abs(X_numerical-X_true)];
   difference = [difference, abs(solution01(t) - X_true)];
end

figure()
% Log-log of local truncation error
loglog(h_list, local_error, '.b', 'MarkerSize', 12)
title("Local Truncation Error, Midpoint")
xlabel("Step Size")
ylabel("Error")
[p,k] = loglog_fit(h_list,local_error, filterparams)
hold on
% Difference Line
loglog(h_list, difference, '.r', 'MarkerSize', 12)
% Fit Line
loglog(h_list, k*h_list.^p, '--m', 'LineWidth', 1.5)
legend("Local Truncation Error", "Difference", "Fit Line")
%% rate_func01
function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
X = cos(t);
end
%fits the relationship y = k*x^p to (x,y) data
%INPUTS:
%x_regression: the x data (a row or column vector)
%y_regression: the y data (a row or column vector)
%filter_params (optional): a struct with parameters to set the range of
%data points that should be included in the regression
%you can set as many or as few of these parameters as you'd like
%   filter_params.min_index
%   filter_params.max_index
%   filter_params.min_xval
%   filter_params.max_xval
%   filter_params.min_yval
%   filter_params.max_yval
%OUTPUTS:
%[p,k]: the regressed values for relationship y = k*x^p
function [p,k] = loglog_fit(x_regression,y_regression,varargin)
    
    %convert x_regression to a column vector if it's a row vector
    if size(x_regression,1)==1
        x_regression = abs(x_regression)';
    end

    %convert y_regression to a column vector if it's a row vector
    if size(y_regression,1)==1
        y_regression = abs(y_regression)';
    end

    %if filter_params has been provided, then filter the data points
    if nargin==3
        filter_params = varargin{1};
        
        num_points = length(x_regression);
        indices = 1:num_points;

        filter_bool = ones(num_points,1);

        if isfield(filter_params,'min_index')
            filter_bool = filter_bool & indices>=filter_params.min_index;
        end

        if isfield(filter_params,'max_index')
            filter_bool = filter_bool & indices<=filter_params.max_index;
        end

        if isfield(filter_params,'min_xval')
            filter_bool = filter_bool & x_regression>=filter_params.min_xval;
        end

        if isfield(filter_params,'max_xval')
            filter_bool = filter_bool & x_regression<=filter_params.max_xval;
        end

        if isfield(filter_params,'min_yval')
            filter_bool = filter_bool & y_regression>=filter_params.min_yval;
        end

        if isfield(filter_params,'max_yval')
            filter_bool = filter_bool & y_regression<=filter_params.max_yval;
        end

        x_regression = x_regression(filter_bool);
        y_regression = y_regression(filter_bool);
    end

    %compute the logs of x_regression and y_regression
    Y = log(y_regression);
    X1 = log(x_regression);

    %set up the regression
    X2 = ones(length(X1),1);

    %run the regression
    coeff_vec = regress(Y,[X1,X2]);
    
    %pull out the coefficients from the fit
    p = coeff_vec(1);
    k = exp(coeff_vec(2));
end