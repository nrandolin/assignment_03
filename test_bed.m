clear

%% Test Forward Euler
tspan = [0,5];
X0 = [1; 0];
h_ref = 0.01;
[t_list_01,X_list_01,h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);

t = 0:0.001:5;
X_true = solution01(t);

figure()
plot(t, X_true)
hold on
plot(t_list_01, X_list_01)
title("Forward Euler (h=0.01)")
legend("True Solution", "Forward Euler")
xlabel("t")
ylabel("X")
hold off

tspan = [0,5];
X0 = [1; 0];
h_ref = 0.1;

[t_list_01,X_list_01,h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
figure()
plot(t, X_true)
hold on
plot(t_list_01, X_list_01)
title("Forward Euler (h=0.1)")
legend("True Solution", "Forward Euler")
xlabel("t")
ylabel("X")

%% rate_func01
function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
X = cos(t);
end