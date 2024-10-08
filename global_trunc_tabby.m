
%% Test Global h/p val
t0 = 0;          % Start time
tf = 0.5;        % End time (replace with the desired final time)
tspan = [t0, tf];  % Full integration interval
X0 = [0; 1];      % Initial conditions

h_list = logspace(-5,1,50);  % Time step sizes
global_error_euler = [];
global_error_mid = [];
rate_function_calls = [];
X_true = [];

for i = 1:length(h_list)
   h_ref = h_list(i);
   [t_list, X_list, h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01, tspan, X0, h_ref);
[t_list_mid, X_list_mid, h_avg_mid, num_evals_mid] = explicit_midpoint_fixed_step_integration(@rate_func01, tspan, X0, h_ref);
  
X_numerical = X_list(end); 
X_numerical_mid = X_list_mid(end); 
   X_true = [X_true, solution01(tf)];  

   global_error_euler = [global_error_euler, norm(X_numerical - X_true(end))];
   global_error_mid = [global_error_mid, norm(X_numerical_mid - X_true(end))];
   rate_function_calls = [rate_function_calls, num_evals];
end

[p,k] = loglog_fit(h_list,global_error)

% Plotting Global Truncation Error vs Step Size
figure;
loglog(h_list, global_error_euler, 'o-'); hold on;
xlabel('Step size (h)');
ylabel('Global Truncation Error');
title('Global Truncation Error vs Step Size');
grid on; hold on;
loglog(h_list, global_error_mid, 'o-', 'Color','r');


hold off;
%% Plotting Global Error vs Rate Function Calls
figure;
loglog(rate_function_calls, global_error_euler, 'o-');
xlabel('Number of rate function calls');
ylabel('Global Truncation Error');
title('Global Error vs Number of Rate Function Calls Forward Euler Algorithm');
grid on;hold on;
loglog(rate_function_calls, global_error_mid, 'o-', 'Color','r');
legend()

%% rate_func01
function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
X = cos(t);
end

