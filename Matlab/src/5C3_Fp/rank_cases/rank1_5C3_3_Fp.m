function [u_sol,v_sol] = rank1_5C3_3_Fp(r,prime)
%RANK1_5C3_3_FP Solve the resulting subsystem of equations for the case R3.
arguments
    r (1,4) {mustBeReal}
    prime (1,1) {mustBePrime}
end
% Obtain results from equations
v_sol = (0:prime-1)';
u_sol = FF(-r(1,3).*v_sol-r(1,4),prime).value;
end