function [v_sol,w_sol] = rank4_5C3_1234_Fp(r,prime)
%RANK4_5C3_1234_FP Solve the resulting subsystem of equations for the case R1234.
arguments
    r (4,1) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePositive,mustBeInteger}
end
% Obtain solutions from equations
v_root = FF(-r(3,1),prime).value;
w_root = FF(-r(4,1),prime).value;
% Remove roots not satysfying the conditions
cond = (FF(-r(1,1),prime) == FF(r(3,1),prime)^2) & ...
    (FF(-r(2,1),prime) == FF(r(4,1),prime)^2);
v_sol = v_root(cond);
w_sol = w_root(cond);
end