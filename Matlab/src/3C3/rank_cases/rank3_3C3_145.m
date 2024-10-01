function [v_sol,w_sol] = rank3_3C3_145(r)
%RANK_3C3_145 Solve the resulting subsystem of equations for the case R145.
arguments
    r (3,3) {mustBeReal}
end
% Obtain root candidates from equations
w_root = -r(3,3);
v_root = -r(2,3);
% Remove roots not satysfying the conditions
conds = abs(-r(1,1).*v_root.^2-r(1,2).*w_root.^2-r(1,3)-w_root.^3) < 10^-10;
v_sol = v_root(conds);
w_sol = w_root(conds);
end