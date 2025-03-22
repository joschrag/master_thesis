function [v_sol,w_sol] = rank3_3C3_245(r)
%RANK_3C3_245 Solve the resulting subsystem of equations for the case R245.
arguments
    r (3,3) 
end
% Obtain root candidates from equations
v_root = -r(2,3);
w_root = -r(3,3);
% Remove roots not satysfying the conditions
conds = abs(-r(1,2).*w_root.^2-r(1,3)-v_root.^2) < 10^-10;
v_sol = v_root(conds);
w_sol = w_root(conds);
end