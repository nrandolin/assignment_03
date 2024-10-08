clear all
t = 0.5;
XA = solution01(t);
h=0.1;
h_list = linspace(10E-5,10,50);
tspan = [0,t];
X0 = [0; 1];
h_ref = 0.1;
[XB,num_evals] = implicit_midpoint_step(@rate_func01,t,XA,h);
[XB,num_evals] = backward_euler(@rate_func01,t,XA,h);
[t_list,X_list,h_avg, num_evals] = fixed_step_integration(@rate_func01,@forward_euler_step,tspan,X0,h_ref);

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
 
    X_list = zeros(num_steps+1,length(X0));
    X_list(1,:) = X0';


    for i = 1:num_steps
        t = t_list(i);
        [XB, temp_eval] = step_func(rate_func_in,t,XA,h_avg);
        num_evals = num_evals + temp_eval;

        X_list(i+1,:)= XB';
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

