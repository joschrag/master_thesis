function [v_sol,w_sol] = rank2_5C3_13_Fp(r,prime)
%RANK2_5C3_13_FP Solve the resulting subsystem of equations for the case R13.
arguments
    r (2,3) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePrimeOrZero,mustBePositive}
end
v = sym("v",["real","integer"]);
w = sym("w",["real","integer"]);
% Obtain results from equations
v0 = FF(-r(2,2)*w-r(2,3),prime);
w_pol = -r(1,1)*w^2-r(1,2)*w-r(1,3)-v^2;
w_pol_s = FF(subs(w_pol,v,v0.value),prime);
w_sol = get_gf_root(coeffs(w_pol_s.value,w,"All"),prime);
if ~isempty(w_sol)
    v_sol = FF(-(r(2,2).*w_sol+r(2,3)),prime).value;
else
    v_sol = [];
end
end