function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
X = cos(t);
end

function dXdt = rate_func02(t,X)
dXdt = [0,-1;1,0]*X;
end

function X = solution02(t)
X = [cos(t);sin(t)];
end
