function [v_sol,w_sol] = rank2_5C3_14_Fp(r,prime)
%RANK2_5C3_14_FP Solve the resulting subsystem of equations for the case R14.
arguments
    r (2,3) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePrimeOrZero,mustBePositive}
end
% Obtain solutions from equations
w_root = FF(-r(2,3),prime).value;
v_sol = get_gf_root([-1,-r(1,2),-r(1,1)*r(2,3)^2-r(1,3)],prime);
w_sol = repmat(w_root,size(v_sol));
end