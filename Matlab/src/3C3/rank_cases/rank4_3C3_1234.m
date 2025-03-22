function [v_sol,w_sol] = rank4_3C3_1234(r,complex)
%RANK4_3C3_1234 Solve the subsystem of equations for the case R1234.
arguments
    r (4,2) 
    complex (1,1) {mustBeNumericOrLogical} = false;
end
% Obtain root candidates from equations
w_root = roots([-1,-r(3,1),-r(3,2)]);
if ~complex
w_root = w_root(abs(imag(w_root))<10^-10);
end
v_root = -r(4,1).*w_root-r(4,2);
% Remove roots not satysfying the conditions
conds = abs(-r(1,1).*w_root-r(1,2)-w_root.^3)<10^-10 & abs(-r(2,1).*w_root-r(2,2)-v_root.^2)<10^-10;
v_sol = v_root(conds);
w_sol = w_root(conds);
end

