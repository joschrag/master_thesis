function [u_sol,v_sol] = rank3_5C3_234_Fp(r,prime)
%RANK3_5C3_234_FP Solve the resulting subsystem of equations for the case R234.
arguments
    r (3,2) {mustBeReal}
    prime (1,1) {mustBePrime}
end
% Obtain solutions from equations
u_root = FF(-r(2,2),prime).value;
v_root = FF(-r(3,2),prime).value;
% Remove roots not satysfying the conditions
cond =  FF(-r(1,2)-r(3,2)^2,prime).value == 0;
v_sol = v_root(cond);
u_sol = u_root(cond);
end