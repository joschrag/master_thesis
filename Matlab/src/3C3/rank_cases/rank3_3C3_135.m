function [v_sol,w_sol] = rank3_3C3_135(r)
%RANK_3C3_135 Solve the resulting subsystem of equations for the case R135.
arguments
    r (3,3) {mustBeReal}
end
% Obtain root candidates from equations
w_root = -r(3,3);
if r(2,2) == 0
    v_root = roots([-r(1,1),-r(1,2),-r(1,3)-w_root^3]);
    v_root = v_root(abs(imag(v_root))<10^-10);
    if ~isempty(v_root)
        conds = abs(-r(2,3)-w_root.^2) < 10^-10;
        w_root = repmat(w_root,size(v_root));
    else
        conds = [];
    end
else
    v_root = -(w_root.^2+r(2,3))/r(2,2);
    conds = abs(-r(1,1).*v_root.^2-r(1,2).*v_root-r(1,3)-w_root.^3) < 10^-10;
end
% Remove roots not satysfying the conditions
v_sol = v_root(conds);
w_sol = w_root(conds);
end