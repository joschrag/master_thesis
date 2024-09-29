function [v_sol,w_sol] = rank5_4C3_12345(r)
%RANK5_4C3_12345 Solve the resulting subsystem of equations for the case R12345.
arguments
    r (5,1) {mustBeReal}
end
% Obtain solution candidates from equations
v_root = -r(4,1);
w_root = -r(5,1);
% Remove roots not satysfying the conditions
cond =  abs(r(1,1) + r(4,1)^2) < 10^-10 &...
    abs(r(3,1) + r(5,1)^2) < 10^-10 &...
    abs(r(2,1) + r(4,1)*r(5,1)) < 10^-10;
v_sol = v_root(cond);
w_sol =  w_root(cond);
end