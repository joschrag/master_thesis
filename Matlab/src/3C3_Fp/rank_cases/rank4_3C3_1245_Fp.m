function [v_sol,w_sol] = rank4_3C3_1245_Fp(r,prime)
%RANK_3C3_1245_FP Solve the resulting subsystem of equations for the case R1245.
arguments
    r (4,2) {mustBeReal}
    prime (1,1) {mustBeInteger,mustBePositive}
end
% Obtain root candidates from equations
v_root = FF(-r(3,2),prime).value;
w_root = FF(-r(4,2),prime).value;
% Remove roots not satysfying the conditions
conds = FF(-r(2,1).*w_root.^2-r(2,2)-v_root.^2,prime).value == 0 & FF(-r(1,1).*w_root.^2-r(1,2)-w_root.^3,prime).value == 0;
v_sol = v_root(conds);
w_sol = w_root(conds);
end

