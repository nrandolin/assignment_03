filterparams.min_xval = 50;
filterparams.max_xval = 1000;

%% Test Global h/p val
t0 = 0;          % Start time
tf = 0.5;        % End time (replace with the desired final time)
tspan = [t0, tf];  % Full integration interval
X0 = 1;      

h_list = logspace(-5,1,100);  % Time step sizes
global_error_euler = [];
global_error_mid = [];
rate_function_calls_mid = [];
rate_function_calls_euler = [];
X_true = [];
results_matrix = zeros(length(h_list), 6);

for i = 1:length(h_list)
   h_ref = h_list(i);
   [t_list, X_list, h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01, tspan, X0, h_ref);
[t_list_mid, X_list_mid, h_avg_mid, num_evals_mid] = explicit_midpoint_fixed_step_integration(@rate_func01, tspan, X0, h_ref);
  
    X_numerical = X_list(end,:)'; 
    X_analytical = solution01(tf);
    X_numerical_mid = X_list_mid(end,:)'; 
   X_true = [X_true, X_analytical];  


   global_error_euler = [global_error_euler, norm(X_numerical - X_analytical)];
   global_error_mid = [global_error_mid, norm(X_numerical_mid - X_analytical)];
   rate_function_calls_euler = [rate_function_calls_euler, num_evals];
   rate_function_calls_mid = [rate_function_calls_mid, num_evals_mid];

end

[p_euler,k_euler] = loglog_fit(rate_function_calls_euler,global_error_euler, filterparams);
[p_mid,k_mid] = loglog_fit(rate_function_calls_mid,global_error_mid, filterparams);

[p_euler_step,k_euler_step] = loglog_fit(h_list,global_error_euler, filterparams);
[p_mid_step,k_mid_step] = loglog_fit(h_list,global_error_mid, filterparams);

    figure;

    % Euler local error data points and fit line
    loglog(rate_function_calls_euler, global_error_euler, 'b', 'LineWidth', 1.5); hold on;
    loglog(rate_function_calls_mid, global_error_mid,'g', 'LineWidth', 1.5, 'DisplayName', 'Midpoint Error');
    
    loglog(rate_function_calls_euler, k_euler*rate_function_calls_euler.^p_euler, 'r--', 'LineWidth', 1.5);
    loglog(rate_function_calls_mid, k_mid*rate_function_calls_mid.^p_mid, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Midpoint Fit Line');
    xlabel('Rate of Function Calls');
    ylabel('Global Truncation Error');
    title('Global Truncation Error vs Rate of Function Calls');
    legend("Euler","Midpoint", "Fit Line")
    hold off;

    figure;

    % Euler local error data points and fit line
    loglog(h_list, global_error_euler, 'b', 'LineWidth', 1.5); hold on;
    loglog(h_list, global_error_mid,'g', 'LineWidth', 1.5, 'DisplayName', 'Midpoint Error');
    
    loglog(h_list, k_euler_step*h_list.^p_euler_step, 'r--', 'LineWidth', 1.5);
    loglog(h_list, k_mid_step*h_list.^p_mid_step, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Midpoint Fit Line');
    xlabel('Step size');
    ylabel('Global Truncation Error');
    title('Global Truncation Error vs Step Size');
    legend("Euler","Midpoint", "Fit Line")
    hold off;


%% rate_func01
function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
X = cos(t);
end

