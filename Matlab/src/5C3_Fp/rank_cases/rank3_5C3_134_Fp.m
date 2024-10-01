function [v_sol,w_sol] = rank3_5C3_134_Fp(r,prime)
%RANK3_5C3_134_FP Solve the resulting subsystem of equations for the case R134.
arguments
    r (3,2) {mustBeReal}
    prime (1,1) {mustBePrimeOrZero,mustBePositive}
end
% Obtain solutions from equations
w_root = FF(-r(3,2),prime).value;
v_root = FF(-r(2,2),prime).value;
% Remove roots not satysfying the conditions
cond = FF(-r(1,1)*r(3,2)^2-r(1,2)-r(2,2)^2,prime).value == 0;
v_sol = v_root(cond);
w_sol = w_root(cond);
end