function [v_sol,w_sol] = rank4_4C3_1234(r)
%RANK4_4C3_1234 Solve the resulting subsystem of equations for the case R1234.
arguments
    r (4,2) {mustBeReal}
end
% Obtain solution candidates from equations
w_root = roots([-1,-r(3,1),-r(3,2)]);
w_root = w_root(abs(imag(w_root))<10^-10);
v_root = -r(4,1).*w_root-r(4,2);
% Remove roots not satysfying the conditions
conds = abs(-r(1,1).*w_root-r(1,2)-v_root.^2) < 10^-10 & abs(-r(2,1).*w_root-r(2,2)-v_root.*w_root) < 10^-10;
v_sol = v_root(conds);
w_sol = w_root(conds);
end