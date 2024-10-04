%% Forward Euler
%Runs numerical integration using forward Euler approximation
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
forward_euler_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
    num_steps = ceil((tspan(2)-tspan(1))/h_ref); % calculate total number of steps
    h_avg = (tspan(2)-tspan(1))/num_steps; % calculate h average
    % define variables
    XA = X0(2);
    num_evals = 0;
    t_list = [];
    t = tspan(1);
    X_list = [];
    %calculate the values until it is just short of the end value
    for i = 1:num_steps-1
        [XB, ~] = forward_euler_step(rate_func_in,t,XA,h_ref);
        num_evals = num_evals + 1;
        t_list = [t_list, t];
        X_list = [X_list, XB];
        XA = XB;
        t = t+h_ref;
    end
    %calculate the final step
    [XB, ~] = forward_euler_step(rate_func_in,t,X_list(end),(tspan(2)-t_list(end)));
    t_list = [t_list, t];
    X_list = [X_list, XB];
end
%% X AT NEXT TIME STEP FORWARD EULER
%This function computes the value of X at the next time step
%using the Forward Euler approximation
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
function [XB,num_evals] = forward_euler_step(rate_func_in,t,XA,h)
    XB = XA + h*rate_func_in(t, XA);
    num_evals = 1;
end
