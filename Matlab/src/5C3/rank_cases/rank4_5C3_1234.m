function [v_sol,w_sol] = rank4_5C3_1234(r)
%RANK4_5C3_1234 Solve the resulting subsystem of equations for the case R1234.
arguments
    r (4,1) {mustBeReal}
end
% Obtain solutions from equations
if abs(r(1,1) + r(3,1)^2) < 10^-10 && ...
        abs(r(2,1) + r(4,1)^2) < 10^-10
    v_sol = -r(3,1);
    w_sol = -r(4,1);
else
    v_sol = [];
    w_sol = [];
end
end