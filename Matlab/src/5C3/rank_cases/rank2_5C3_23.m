function [v_sol,w_sol] = rank2_5C3_23(r)
%RANK2_5C3_23 Solve the resulting subsystem of equations for the case R23.
arguments
    r (2,3) {mustBeReal}
end
%Obtain solutions from equations
w_0 = roots([-1,-r(1,2),-r(1,3)]);
w_0 = w_0(imag(w_0)==0);
if ~isempty(w_0)
    w_sol = w_0;
    v_sol = -r(2,2).*w_0-r(2,3);
else
    w_sol = [];
    v_sol = [];
end
end