function [v_sol,w_sol] = rank4_4C3_1345_Fp(r,prime)
%RANK4_4C3_1345_FP Solve the resulting subsystem of equations for the case R1345.
arguments
    r (4,2) {mustBeReal}
    prime (1,1) {mustBePrimeOrZero}
end
% Obtain solution candidates from equations
v_root = FF(-r(3,2),prime).value;
w_root =  FF(-r(4,2),prime).value;
% Remove roots not satysfying the conditions
conds = FF(-r(1,1).*w_root.*v_root-r(1,2)-v_root.^2,prime).value==0 & FF(-r(2,2)-w_root.^2,prime).value==0;
v_sol = v_root(conds);
w_sol = w_root(conds);
end