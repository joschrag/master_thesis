function [u_sol,v_sol] = rank2_5C3_14_Fp(r,prime)
%RANK2_5C3_14_FP Solve the resulting subsystem of equations for the case R14.
arguments
    r (2,3) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePrime}
end
% Obtain solutions from equations
v_root = FF(-r(2,3),prime).value;
u_sol = get_gf_root([-1,-r(1,2),-r(1,1)*r(2,3)^2-r(1,3)],prime);
v_sol = repmat(v_root,size(u_sol));
end