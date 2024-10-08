%This function computes the value of X at the next time step
%using the implicit midpoint approximation
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
function [XB,num_evals] = implicit_midpoint_step(rate_func_in,t,XA,h)
    func = @(XB) XA + h*rate_func_in(t + (h/2), 0.5*(XA + XB)) - XB;
    [XB, num_evals] = multi_newton_solver(func, XA, true);
end

%%
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
