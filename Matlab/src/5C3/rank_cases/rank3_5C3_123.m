function [v_sol,w_sol] = rank3_5C3_123(r)
%RANK3_5C3_123 Solve the resulting subsystem of equations for the case R123.
arguments
    r (3,2) {mustBeReal}
end
% Obtain solutions from equations
w_sol = roots([-1,-r(2,:)]);
w_sol = unique(w_sol(imag(w_sol)==0));
if isempty(w_sol)
    v_sol = [];
    return
end
v_sol = -(r(3,1).*w_sol+r(3,2));

control = abs(-(r(1,1).*w_sol+v_sol.^2+r(1,2))) < 10^-10;
w_sol = w_sol(control);
v_sol = v_sol(control);
end

