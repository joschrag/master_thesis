function [v_sol,w_sol] = rank3_3C3_125(r)
%RANK_3C3_125 Solve the resulting subsystem of equations for the case R125.
arguments
    r (3,3) {mustBeReal}
end
% Obtain root candidates from equations
w_root = -r(3,3);
v_root = roots([-1,-r(2,2),-r(2,1)*w_root^2-r(2,3)]);
v_root = v_root(abs(imag(v_root))<10^-10);
w_root = repmat(w_root,size(v_root));
% Remove roots not satysfying the conditions
conds = abs(-r(1,1).*w_root.^2-r(1,2).*v_root-r(1,3)-w_root.^3) <10^-10;
v_sol = v_root(conds);
w_sol = w_root(conds);
end