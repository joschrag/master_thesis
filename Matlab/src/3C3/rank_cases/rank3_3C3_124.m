function [v_sol,w_sol] = rank3_3C3_124(r,complex)
%RANK_3C3_124 Solve the resulting subsystem of equations for the case R124.
arguments
    r (3,3) 
    complex (1,1) {mustBeNumericOrLogical} = false;
end
% Obtain root candidates from equations
w_root = roots([-1,-r(1,1),-r(1,2),-r(1,3)]);
if ~complex
    w_root = w_root(abs(imag(w_root))<10^-10);
end
v_root = -r(3,2).*w_root-r(3,3);
% Remove roots not satysfying the conditions
conds = abs(-r(2,1).*w_root.^2-r(2,2).*w_root-r(2,3)-v_root.^2) < 10^-10;
v_sol = v_root(conds);
w_sol = w_root(conds);
end