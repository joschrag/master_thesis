function [v_sol,w_sol] = rank2_5C3_34(r)
%RANK2_5C3_34 Solve the resulting subsystem of equations for the case R34.
arguments
    r (2,3) {mustBeReal}
end
% Obtain solutions from equations
v_sol = -r(1,3);
w_sol = -r(2,3);
end