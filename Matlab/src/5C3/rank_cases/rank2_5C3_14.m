function [v_sol,w_sol] = rank2_5C3_14(r)
%RANK2_5C3_14 Solve the resulting subsystem of equations for the case R14.
arguments
    r (2,3) {mustBeReal}
end
% Obtain solutions from equations
v_sol = [];
w_sol = [];
w_0 = -r(2,3);
v_0 = roots([-1,-r(1,2),-r(1,1)*r(2,3)^2-r(1,3)]);
v_0 = v_0(imag(v_0)==0);
if ~isempty(v_0)
    v_sol = v_0;
    w_sol = repmat(w_0,size(v_0));
end
end