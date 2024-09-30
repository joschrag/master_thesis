function [v_sol,w_sol] = rank1_5C3_4_Fp(r,prime)
%RANK1_5C3_4_FP Solve the resulting subsystem of equations for the case R4.
arguments
    r (1,4) {mustBeReal}
    prime (1,1) {mustBePrime}
end
% Obtain results from equations
v_sol = (0:prime-1)';
w_sol = repmat(FF(-r(1,4),prime).value,prime,1);
end