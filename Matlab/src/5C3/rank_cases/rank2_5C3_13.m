function [v_sol,w_sol] = rank2_5C3_13(r)
%RANK2_5C3_13 Solve the resulting subsystem of equations for the case R13.
arguments
    r (2,3) {mustBeReal}
end
% Obtain solutions from equations
v = sym("v","real");
w = sym("w","real");
v0 = -r(2,2)*w-r(2,3);
w_pol = subs(-r(1,1)*w^2-r(1,2)*w-r(1,3)-v^2,v,v0);
w_sol = roots(coeffs(w_pol,w,"All"));
w_sol = w_sol(imag(w_sol)==0);
if ~isempty(w_sol)
    v_sol = -(r(2,2).*w_sol+r(2,3));
else
    v_sol = [];
end
end