function [v_sol,w_sol] = rank3_5C3_124(r)
%RANK3_5C3_124 Solve the resulting subsystem of equations for the case R124.
arguments
    r (3,2) {mustBeReal}
end
% Obtain solutions from equations
w_0 = -r(3,2);
v_0 = roots([-1,-r(1,1),-r(1,2)]);
v_0 = v_0(imag(v_0)==0);
control = abs(-w_0^2-r(2,1).*v_0-r(2,2)) < 10^-10;
if ~any(control)
    w_sol = [];
    v_sol = [];
    return
end
v_sol = v_0(control);
w_sol = repmat(w_0,size(v_0));
end

