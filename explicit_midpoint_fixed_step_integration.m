%Runs numerical integration using explicit midpoint approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0’;X1’;X2’;...;(X_end)’] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = ...
explicit_midpoint_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
    [num_steps, h_avg] = iteration_solver(tspan, h_ref);
    % define variables
    XA = X0(2);
    num_evals = 0;
    t_list = [];
    t = tspan(1);
    X_list = [];
    %calculate the values until it is just short of the end value
    for i = 1:num_steps
        [XB, ~] = explicit_midpoint_step(rate_func_in,t,XA,h_avg);
        num_evals = num_evals + 1;
        t_list = [t_list, t];
        X_list = [X_list, XB];
        XA = XB;
        t = t+h_ref;
    end
end

%% Single Step
%This function computes the value of X at the next time step
%using the explicit midpoint approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB,num_evals] = explicit_midpoint_step(rate_func_in,t,XA,h)
    X_n_half = XA + h/2*rate_func_in(t, XA);
    XB = X_n_half + h*rate_func_in(t, X_n_half);
    num_evals = 2;
end

%% iteration solver (find actual n and h)
function [num_steps, h] = iteration_solver(tspan, h_ref)
    range = tspan(2)-tspan(1);
    num_steps = range/h_ref;%The number of steps is the range devided by goal h 
    num_steps = ceil(num_steps);%Round the number of steps up (to get a real number)
    h = range/num_steps; % Divide range by steps to get real h
end