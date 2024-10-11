
function local_trunc_tabby()

    %% LOCAL FORWARD EULER
  
h_list = logspace(-5,1,100); 

    local_error_euler = [];
    local_error_mid = [];
    X_true = [];
    difference = [];

    t_ref = 4.49;
    XA = solution01(t_ref);

    for i = 1:length(h_list)
        h_ref = h_list(i);
        
        [XB_euler,~] = forward_euler_step(@rate_func01,t_ref,XA,h_ref);
        [XB_mid,~] = explicit_midpoint_step(@rate_func01,t_ref,XA,h_ref);
 
        X_analytical = solution01(t_ref+h_ref);
        
       X_true = [X_true, X_analytical];
       local_error_euler = [local_error_euler, norm(XB_euler - X_analytical)];
       local_error_mid = [local_error_mid, norm(XB_mid - X_analytical)];
       difference = [difference, norm(X_analytical-solution01(t_ref))];

    end

    
    [p_euler,k_euler] = loglog_fit(h_list,local_error_euler)
    [p_mid,k_mid] = loglog_fit(h_list,local_error_mid)
    [p_diff,k_diff] = loglog_fit(h_list,difference)

    %%%%%
    
    % Plotting Local Truncation Error vs Step Size
    figure; 

    % Euler local error data points and fit line
    loglog(h_list, local_error_euler, 'b', 'LineWidth', 1.5); hold on;
    grid on;
    xlabel('Step size (h)');
    ylabel('Local Truncation Error');

    title('Local Truncation Error vs Step Size');
    loglog(h_list, local_error_mid,'g', 'LineWidth', 1.5, 'DisplayName', 'Midpoint Error');


    loglog(h_list, difference, 'c', 'LineWidth', 1.5, 'DisplayName', 'Difference');
    
    loglog(h_list, k_diff*h_list.^p_diff, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Difference Fit Line');
    loglog(h_list, k_euler*h_list.^p_euler, 'r--', 'LineWidth', 1.5);
        loglog(h_list, k_mid*h_list.^p_mid, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Midpoint Fit Line');

    legend("Euler","Midpoint", "Difference", "Fit Line")
    hold off;
    
    
    %% rate_func01
    function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
    end
    
    function X = solution01(t)
    X = cos(t);
    end
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

function [XB,num_evals] = forward_euler_step(rate_func_in,t,XA,h)
    XB = XA + h*rate_func_in(t, XA);
    num_evals = 1;
end

function [XB,num_evals] = explicit_midpoint_step(rate_func_in,t,XA,h)
    X_n_half = XA + h/2*rate_func_in(t, XA);
    XB = XA + h*rate_func_in(t+h/2, X_n_half);
    num_evals = 2;
end

%% rate_func01
function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
X = cos(t);
end
