%Runs fixed step numerical integration
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%step_func: the function used to approximate X(t) at the next step
% [XB,num_evals] = step_func(rate_func_in,t,XA,h)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0’;X1’;X2’;...;(X_end)’] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = ...
fixed_step_integration_old(rate_func_in,step_func,tspan,X0,h_ref)
    % calculate steps and h
    [num_steps, h_avg] = iteration_solver(tspan, h_ref);
    % define variables
    XA = X0;
    num_evals = 0;
    t_list = linspace(tspan(1),tspan(2),num_steps+1);
 
    X_list = zeros(num_steps+1,length(X0));
    X_list(1,:) = X0';

    % forward euler
    if step_func == 1
        for i = 1:num_steps
            t = t_list(i);
            [XB, temp_eval] = forward_euler_step(rate_func_in,t,XA,h_avg);
            num_evals = num_evals + temp_eval;
    
            X_list(i+1,:)= XB';
            XA = XB;
        end
    end

    % forward midpoint
    if step_func == 2
        for i = 1:num_steps
            t = t_list(i);
            [XB, temp_eval] = explicit_midpoint_step(rate_func_in,t,XA,h_avg);
            num_evals = num_evals + temp_eval;
    
            X_list(i+1,:)= XB';
            XA = XB;
        end
    end

    % Backward Euler
     if step_func == 3
        for i = 1:num_steps
            t = t_list(i);
            [XB, temp_eval] = backward_euler(rate_func_in,t,XA,h_avg);
            num_evals = num_evals + temp_eval;
    
            X_list(i+1,:)= XB';
            XA = XB;
        end
     end
end
%% Testing
function [t_list,X_list,h_avg, num_evals] = ...
fixed_step_integration(rate_func_in,step_func,tspan,X0,h_ref)
    % calculate steps and h
    [num_steps, h_avg] = iteration_solver(tspan, h_ref);
    % define variables
    XA = X0;
    num_evals = 0;
    t_list = linspace(tspan(1),tspan(2),num_steps+1);
 
    X_list = zeros(num_steps+1)
    X_list(1,:) = X0(1)


    for i = 1:num_steps
        t = t_list(i);
        [XB, temp_eval] = step_func(rate_func_in,t,XA,h_avg);
        num_evals = num_evals + temp_eval;

        X_list(i+1,:)= XB(1);
        XA = XB;
    end
end


%% Single Step midpoint
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
    XB = XA + h*rate_func_in(t+h/2, X_n_half);
    num_evals = 2;
end
%% Single Step EULER
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

%% iteration solver (find actual n and h)
function [num_steps, h] = iteration_solver(tspan, h_ref)
    range = tspan(2)-tspan(1);
    num_steps = range/h_ref;%The number of steps is the range devided by goal h 
    num_steps = ceil(num_steps);%Round the number of steps up (to get a real number)
    h = range/num_steps; % Divide range by steps to get real h
end
%This function computes the value of X at the next time step
%using the Backward Euler approximation
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
function [XB,num_evals] = backward_euler(rate_func_in,t,XA,h)
    func = @(XB) XA + h*rate_func_in(t+h, XB) - XB;
    [XB, num_evals] = multi_newton_solver(func, XA, true);
end

%% MULTI NEWTON SOLVER
function [x_next, num_evals] = multi_newton_solver(fun,x_guess,varargin)
%true if supposed to use analytical jacobian, false otherwise
use_analytical_jacobian = nargin==3 && varargin{1}(1);

A_thresh = 10e-14;
B_thresh = 10e-14;
num_evals = 0;

    if use_analytical_jacobian == true
    f_val = fun(x_guess);
    num_evals = 1;
    J = approximate_jacobian(fun,x_guess);
    x_next = x_guess;
        while max(abs(J\f_val))> A_thresh && max(abs(f_val))>B_thresh && abs(det(J*J')) > 1e-14
            %calculate next x
            x_next = x_next - J\f_val;
            %establish next y value
            f_val = fun(x_next);
            J = approximate_jacobian(fun,x_next);  
            num_evals = num_evals+1;
        end
    end

%Loop through until x is small
    if use_analytical_jacobian == false
    [f_val, J] = fun(x_guess);
    x_next = x_guess;
        while max(abs(J\f_val))> A_thresh && max(abs(f_val))>B_thresh && abs(det(J*J')) > 1e-14
            %calculate next x
            x_next = x_next - J\f_val;
            %establish next y value
            [f_val, J]  = fun(x_next);        
        end
    end

end
%% APPROX JACOBIAN
function J = approximate_jacobian(fun,x)
    % Set initial variables
    ej = zeros(length(x), 1); %variable to store vector of multiplyers
    h = 1e-6;
    J = zeros(length(fun(x)), length(x));
    for i = 1:size(J, 2)
        ej(i) = 1;
        % calculate the partial derivative 
        J(:, i) = (fun(x+h*ej) - fun(x-h*ej))/(2*h);
        ej(i) = 0;
    end 
end