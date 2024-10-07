%% Test iteration solver
tspan = [0,4.9];
h_ref = 0.2;
[num_steps, h] = iteration_solver(tspan, h_ref)
t_final = num_steps*h

%% Test Forward Euler
clear all
tspan = [0,4.9];
X0 = [0; 1];
h_ref = 0.2;
[t_list_01,X_list_01,h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);

t = 0:0.001:5;
X_true = solution01(t);

figure()
plot(t, X_true)
hold on
plot(t_list_01, X_list_01)
title("Forward Euler")
xlabel("t")
ylabel("X")

tspan = [0,5];
X0 = [1; 0];
h_ref = 0.1;

[t_list_02,X_list_02,h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
plot(t_list_02, X_list_02)
legend("True Solution", "h = 0.2", "h = 0.1")
hold off


%% test midpoint
clear all

tspan = [0,4.9];
X0 = [0; 1];
h_ref = 0.2;

[t_list_01,X_list_01,h_avg, num_evals] = explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);

t = 0:0.001:5;
X_true = solution01(t);

figure()
plot(t, X_true)
hold on
plot(t_list_01, X_list_01)
title("Midpoint (h=0.01)")
xlabel("t")
ylabel("X")


tspan = [0,4.9];
X0 = [1; 0];
h_ref = 0.1;

[t_list_02,X_list_02,h_avg, num_evals] = explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);

plot(t_list_02, X_list_02)
title("Midpoint")
legend("True Solution", "h=0.2", "h=0.1")
xlabel("t")
ylabel("X")

%% rate_func01
function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
X = cos(t);
end

%% iteration solver (find actual n and h)
function [num_steps, h] = iteration_solver(tspan, h_ref)
    range = tspan(2)-tspan(1);
    num_steps = range/h_ref;%The number of steps is the range devided by goal h 
    num_steps = ceil(num_steps);%Round the number of steps up (to get a real number)
    h = range/num_steps; % Divide range by steps to get real h
end