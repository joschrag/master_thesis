function [v_sol,w_sol] = rank5_3C3_12345_Fp(r,prime)
%RANK5_3C3_12345_FP Solve the resulting subsystem of equations for the case R12345.
arguments
    r (5,1) {mustBeReal}
    prime (1,1) {mustBePrime}
end
% Obtain root candidates from equations
v_root = FF(-r(4,1),prime).value;
w_root = FF(-r(5,1),prime).value;
% Remove roots not satysfying the conditions
conds =  FF(r(1,1) - r(5,1)^3,prime).value == 0 &...
    FF(r(3,1) + r(5,1)^2,prime).value == 0 &...
    FF(r(2,1) + r(4,1)^2,prime).value == 0;
v_sol = v_root(conds);
w_sol = w_root(conds);
end

