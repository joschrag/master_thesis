function [v_sol,w_sol] = rank4_3C3_1245(r)
%RANK4_3C3_1245 Solve the resulting subsystem of equations for the case R1245.
arguments
    r (4,2) 
end
% Obtain root candidates from equations
v_root = -r(3,2);
w_root = -r(4,2);
% Remove roots not satysfying the conditions
conds = abs(-r(2,1).*w_root.^2-r(2,2)-v_root.^2)< 10^-10 & abs(-r(1,1).*w_root.^2-r(1,2)-w_root.^3)< 10^-10;
v_sol = v_root(conds);
w_sol = w_root(conds);
end

