function [v_sol,w_sol] = rank4_4C3_2345(r)
%RANK4_4C3_2345 Solve the resulting subsystem of equations for the case R123.
arguments
    r (4,2) {mustBeReal}
end
% Obtain solution candidates from equations
v_root = -r(3,2);
w_root = -r(4,2);
% Remove roots not satysfying the conditions
conds = abs(-r(1,2)-w_root.*v_root) < 10^-10 & abs(-r(2,2)-w_root.^2) < 10^-10;
v_sol = v_root(conds);
w_sol = w_root(conds);
end