function [u_sol,v_sol] = rank2_5C3_34_Fp(r,prime)
%RANK2_5C3_34_FP Solve the resulting subsystem of equations for the case R34.
arguments
    r (2,3) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePrime}
end
% Obtain solutions from equations
u_sol = FF(-r(1,3),prime).value;
v_sol = FF(-r(2,3),prime).value;
end