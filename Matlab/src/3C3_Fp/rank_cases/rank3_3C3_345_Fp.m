function [v_sol,w_sol] = rank3_3C3_345_Fp(r,prime)
%RANK_3C3_345_FP Solve the resulting subsystem of equations for the case R345.
arguments
    r (3,3) {mustBeReal}
    prime (1,1) {mustBePrime}
end
% Obtain root candidates from equations
v_root = FF(-r(2,3),prime).value;
w_root = FF(-r(3,3),prime).value;
% Remove roots not satysfying the conditions
conds = FF(-r(1,3)-w_root.^2,prime).value ==0;
v_sol = v_root(conds);
w_sol = w_root(conds);
end

