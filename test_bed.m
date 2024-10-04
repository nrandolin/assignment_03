%% Test Forward Euler
tspan = [0,5];
X0 = [1; 0];
h_ref = 0.01;

[t_list,X_list,h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);


%% rate_func01
function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end