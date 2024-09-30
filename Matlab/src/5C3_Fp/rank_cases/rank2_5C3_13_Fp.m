function [u_sol,v_sol] = rank2_5C3_13_Fp(r,prime)
%RANK2_5C3_13_FP Solve the resulting subsystem of equations for the case R13.
arguments
    r (2,3) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePrime}
end
u = sym("t",["real","integer"]);
v = sym("u",["real","integer"]);
% Obtain results from equations
t0 = FF(-r(2,2)*v-r(2,3),prime);
v_pol = -r(1,1)*v^2-r(1,2)*v-r(1,3)-u^2;
v_pol_s = FF(subs(v_pol,u,t0.value),prime);
v_sol = get_gf_root(coeffs(v_pol_s.value,v,"All"),prime);
if ~isempty(v_sol)
    u_sol = FF(-(r(2,2).*v_sol+r(2,3)),prime).value;
else
    u_sol = [];
end
end