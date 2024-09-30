function [v_sol,w_sol] = rank3_5C3_234(r)
%RANK3_5C3_234 Solve the resulting subsystem of equations for the case R234.
arguments
    r (3,2) {mustBeReal}
end
% Obtain solutions from equations
if abs(-r(1,2)-r(3,2)^2) < 10^-10
    w_sol = -r(3,2);
    v_sol = -r(2,2);
else
    w_sol = [];
    v_sol = [];
    return
end
end

