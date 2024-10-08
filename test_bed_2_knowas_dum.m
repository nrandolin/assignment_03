clear all
t = 0.5;
XA = solution01(t);
h=0.1;

[XB,num_evals] = backward_euler(@rate_func01,t,XA,h);

%% rate_func01
function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
X = cos(t);
end
