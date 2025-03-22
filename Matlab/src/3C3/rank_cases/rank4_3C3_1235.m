function [v_sol,w_sol] = rank4_3C3_1235(r,complex)
%RANK4_3C3_1235 Solve the resulting subsystem of equations for the case R1235.
arguments
    r (4,2) 
    complex (1,1) {mustBeNumericOrLogical} = false;
end
% Obtain root candidates from equations
v_root = roots([-1,-r(2,1),-r(2,2)]);
if ~complex
v_root = v_root(abs(imag(v_root))<10^-10);
end
w_root = repmat(-r(4,2),size(v_root));
% Remove roots not satysfying the conditions
conds = abs(-r(3,1).*v_root-r(3,2)-w_root.^2)<10^-10 & abs(-r(1,1).*v_root-r(1,2)-w_root.^3)<10^-10;
v_sol = v_root(conds);
w_sol = w_root(conds);
end