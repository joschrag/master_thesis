function [v_sol,w_sol] = rank3_3C3_235(r,complex)
%RANK_3C3_235 Solve the resulting subsystem of equations for the case R235.
arguments
    r (3,3) 
    complex (1,1) {mustBeNumericOrLogical} = false;
end
% Obtain root candidates from equations
v_root = roots([-1,-r(1,2),-r(1,3)]);
if ~complex
v_root = v_root(abs(imag(v_root))<10^-10);
end
w_root = repmat(-r(3,3),size(v_root));
% Remove roots not satysfying the conditions
conds = abs(-r(2,2).*v_root-r(2,3)-w_root.^2) <10^-10;
v_sol = v_root(conds);
w_sol = w_root(conds);
end