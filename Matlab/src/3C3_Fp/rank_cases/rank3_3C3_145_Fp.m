function [v_sol,w_sol] = rank3_3C3_145_Fp(r,prime)
%RANK_3C3_145_FP Solve the resulting subsystem of equations for the case R145.
arguments
    r (3,3) {mustBeReal}
    prime (1,1) {mustBeInteger,mustBePositive}
end
% Obtain root candidates from equations
w_root = FF(-r(3,3),prime).value;
v_root = FF(-r(2,3),prime).value;
% Remove roots not satysfying the conditions
conds = FF(-r(1,1).*v_root.^2-r(1,2).*w_root.^2-r(1,3)-w_root.^3,prime).value ==0;
v_sol = v_root(conds);
w_sol = w_root(conds);
end

