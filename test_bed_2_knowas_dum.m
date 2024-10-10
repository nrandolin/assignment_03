%%
clear all
filterparams.min_xval = 0;
filterparams.max_xval = 0.1;

%%
clear all
t = 20;
XA = solution01(t);
h=0.1;
h_list = linspace(10E-5,10,50);
tspan = [0,t];
X0 = [0; 1];
h_ref = 0.38;
% [XB,num_evals] = implicit_midpoint_step(@rate_func01,t,XA,h);
% [XB,num_evals] = backward_euler(@rate_func01,t,XA,h);

subplot(4,1,1);
[t_list2,X_list2,h_avg, num_evals] = fixed_step_integration(@rate_func01,@forward_euler_step,tspan,X0,h_ref)
num_sol = plot(t_list2, X_list2, 'b', 'LineWidth', 1.5, 'DisplayName', 'Numerical Solution')
title("Forward Euler Numerical Solution")
xlabel("Time")
ylabel("X")
hold on
X_true = solution01(t_list2);
true_sol = plot(t_list2, X_true, '--r', 'LineWidth', 1.5, 'DisplayName', 'True Solution')
legend([num_sol(1), true_sol(1)], 'Numerical Solution', 'True Solution')

subplot(4,1,2);
[t_list1,X_list1,h_avg, num_evals] = fixed_step_integration(@rate_func01,@explicit_midpoint_step,tspan,X0,h_ref)
num_sol = plot(t_list1, X_list1, 'b', 'LineWidth', 1.5, 'DisplayName', 'Numerical Solution')
title("Explicit Midpoint Numerical Solution")
xlabel("Time")
ylabel("X")
hold on
X_true = solution01(t_list1);
true_sol = plot(t_list1, X_true, '--r', 'LineWidth', 1.5, 'DisplayName', 'True Solution')
legend([num_sol(1), true_sol(1)], 'Numerical Solution', 'True Solution')

subplot(4,1,3);
[t_list3,X_list3,h_avg, num_evals] = fixed_step_integration(@rate_func01,@implicit_midpoint_step,tspan,X0,h_ref)
num_sol = plot(t_list3, X_list3, 'b', 'LineWidth', 1.5, 'DisplayName', 'Numerical Solution')
title("Implicit Midpoint Numerical Solution")
xlabel("Time")
ylabel("X")
hold on
X_true = solution01(t_list3);
true_sol = plot(t_list3, X_true, '--r', 'LineWidth', 1.5, 'DisplayName', 'True Solution')
legend([num_sol(1), true_sol(1)], 'Numerical Solution', 'True Solution')

subplot(4,1,4);
[t_list4,X_list4,h_avg, num_evals] = fixed_step_integration(@rate_func01,@backward_euler,tspan,X0,h_ref)
num_sol = plot(t_list4, X_list4, 'b', 'LineWidth', 1.5, 'DisplayName', 'Numerical Solution')
title("Backward Euler Numerical Solution")
xlabel("Time")
ylabel("X")
hold on
X_true = solution01(t_list4);
true_sol = plot(t_list4, X_true, '--r', 'LineWidth', 1.5, 'DisplayName', 'True Solution')
legend([num_sol(1), true_sol(1)], 'Numerical Solution', 'True Solution')
%% GLOBAL
t = 0.5;
h_list = linspace(10E-5,10,100);
h_list = logspace(-5,1,50);
tspan = [0,t];
X0 = [0; 1];
glob_error_euler = [];
glob_error_mid = [];
glob_error_backward = [];
glob_error_imp_mid = [];
difference = [];

for i = 1:length(h_list)
   h_ref = h_list(i);
   % Calculate numperical x value
   [t_list1,X_list1,h_avg1, num_evals1] = fixed_step_integration(@rate_func01,@forward_euler_step,tspan,X0,h_ref);
   [t_list2,X_list2,h_avg2, num_evals2] = fixed_step_integration(@rate_func01,@explicit_midpoint_step,tspan,X0,h_ref);
   [t_list3,X_list3,h_avg3, num_evals3] = fixed_step_integration(@rate_func01,@backward_euler,tspan,X0,h_ref);
   [t_list4,X_list4,h_avg4, num_evals4] = fixed_step_integration(@rate_func01,@implicit_midpoint_step,tspan,X0,h_ref);

   X_numerical1 = X_list1(end);
   X_numerical2 = X_list2(end);
   X_numerical3 = X_list3(end);
   X_numerical4 = X_list4(end);
   % calculate the real x value
   X_true = solution01(t+h_ref);

   glob_error_euler = [glob_error_euler, abs(X_numerical1-X_true)];
   glob_error_mid = [glob_error_mid, abs(X_numerical2-X_true)];
   glob_error_backward = [glob_error_backward, abs(X_numerical3-X_true)];
   glob_error_imp_mid = [glob_error_imp_mid, abs(X_numerical4-X_true)];
end
[p_g_euler,k_g_euler] = loglog_fit(h_list,local_error_euler, filterparams)
[p_g_mid,k_g_mid] = loglog_fit(h_list,local_error_mid, filterparams)
[p_g_backward,k_g_backward] = loglog_fit(h_list,local_error_backward, filterparams)
[p_g_imp_mid,k_g_imp_mid] = loglog_fit(h_list,local_error_imp_mid, filterparams)

figure()
% Log-log of local truncation error
loglog(h_list, glob_error_euler, '.b', 'MarkerSize', 10)
hold on
loglog(h_list, glob_error_mid, '.r', 'MarkerSize', 10)
hold on
loglog(h_list, glob_error_backward, '.g', 'MarkerSize', 10)
hold on
loglog(h_list, glob_error_imp_mid, '.m', 'MarkerSize', 10)
title("Local Truncation Error, All Methods")
xlabel("Step Size")
ylabel("Error")
legend("Forward Euler", "Explicit Midpoint", "Backward Euler", "Implicit Midpoint")

%% LOCAL
t0 = 0;          % Start time
tf = 0.5;        % End time (replace with the desired final time)
tspan = [t0, tf];  % Full integration interval
X0 = 1;      % Initial conditions

h_list = logspace(-5,1,50);  % Time step sizes
local_error_euler = [];
local_error_mid = [];
local_error_backward = [];
local_error_imp_mid = [];
rate_function_calls = [];
X_true = [];

t_ref = 4.49;
XA = solution01(t_ref);

for i = 1:length(h_list)
    h_ref = h_list(i);
    
    [XB_euler,~] = forward_euler_step(@rate_func01,t_ref,XA,h_ref);
    [XB_mid,~] = explicit_midpoint_step(@rate_func01,t_ref,XA,h_ref);
    [XB_backward,~] = backward_euler(@rate_func01,t_ref,XA,h_ref);
    [XB_imp_mid,~] = implicit_midpoint_step(@rate_func01,t_ref,XA,h_ref);

    X_analytical = solution01(t_ref+h_ref);
    
    X_true = [X_true, X_analytical];  

   local_error_euler = [local_error_euler, norm(XB_euler - X_analytical)];
   local_error_mid = [local_error_mid, norm(XB_mid - X_analytical)];
   local_error_backward = [local_error_backward, norm(XB_backward - X_analytical)];
   local_error_imp_mid = [local_error_imp_mid, norm(XB_imp_mid - X_analytical)];
   
end
[p_euler,k_euler] = loglog_fit(h_list,local_error_euler, filterparams)
[p_mid,k_mid] = loglog_fit(h_list,local_error_mid, filterparams)
[p_backward,k_backward] = loglog_fit(h_list,local_error_backward, filterparams)
[p_imp_mid,k_imp_mid] = loglog_fit(h_list,local_error_imp_mid, filterparams)

figure()
% Log-log of local truncation error
loglog(h_list, local_error_euler, '.b', 'MarkerSize', 10)
hold on
loglog(h_list, local_error_mid, '.r', 'MarkerSize', 10)
hold on
loglog(h_list, local_error_backward, '.g', 'MarkerSize', 10)
hold on
loglog(h_list, local_error_imp_mid, '.m', 'MarkerSize', 10)
title("Local Truncation Error, All Methods")
xlabel("Step Size")
ylabel("Error")
legend("Forward Euler", "Explicit Midpoint", "Backward Euler", "Implicit Midpoint")
%% rate_func01
function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
X = cos(t);
end

function [XB,num_evals] = forward_euler_step(rate_func_in,t,XA,h)
    XB = XA + h*rate_func_in(t, XA);
    num_evals = 1;
end

function [t_list,X_list,h_avg, num_evals] = ...
fixed_step_integration(rate_func_in,step_func,tspan,X0,h_ref)
    % calculate steps and h
    [num_steps, h_avg] = iteration_solver(tspan, h_ref);
    % define variables
    XA = X0;
    num_evals = 0;
    t_list = linspace(tspan(1),tspan(2),num_steps+1);
 
    X_list = zeros(num_steps+1);
    X_list(1,:) = X0(1);


    for i = 1:num_steps
        t = t_list(i);
        [XB, temp_eval] = step_func(rate_func_in,t,XA,h_avg);
        num_evals = num_evals + temp_eval;

        X_list(i+1,:)= XB(1);
        XA = XB;
    end
end

function [num_steps, h] = iteration_solver(tspan, h_ref)
    range = tspan(2)-tspan(1);
    num_steps = range/h_ref;%The number of steps is the range devided by goal h 
    num_steps = ceil(num_steps);%Round the number of steps up (to get a real number)
    h = range/num_steps; % Divide range by steps to get real h
end

function [XB,num_evals] = backward_euler(rate_func_in,t,XA,h)
    func = @(XB) XA + h*rate_func_in(t+h, XB) - XB;
    [XB, num_evals] = multi_newton_solver(func, XA, true);
end

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
%%
function [XB,num_evals] = explicit_midpoint_step(rate_func_in,t,XA,h)
    X_n_half = XA + h/2*rate_func_in(t, XA);
    XB = XA + h*rate_func_in(t+h/2, X_n_half);
    num_evals = 2;
end


