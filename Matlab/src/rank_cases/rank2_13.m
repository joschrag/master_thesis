function [u_sol,v_sol] = rank2_13(r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (2,3) {mustBeReal}
end
u = sym("u","real");
v = sym("v","real");
t0 = -r(2,2)*v-r(2,3);
v_pol = subs(-r(1,1)*v^2-r(1,2)*v-r(1,3)-u^2,u,t0);
v_sol = roots(coeffs(v_pol,v,"All"));
v_sol = v_sol(imag(v_sol)==0);
if ~isempty(v_sol)
    u_sol = -(r(2,2).*v_sol+r(2,3));
else
    u_sol = [];
end
end