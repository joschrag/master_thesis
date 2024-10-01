function [v_sol,w_sol] = rank3_3C3_135_Fp(r,prime)
%RANK_3C3_135_FP Solve the resulting subsystem of equations for the case R135.
arguments
    r (3,3) {mustBeReal}
    prime (1,1) {mustBePrime}
end
% Obtain root candidates from equations
w_root = FF(-r(3,3),prime).value;
if r(2,2) == 0
    v_root = get_gf_root([-r(1,1),-r(1,2),-r(1,3)-w_root^3],prime);
    if ~isempty(v_root)
        conds = FF(-r(2,3)-w_root.^2,prime).value ==0;
        w_root = repmat(w_root,size(v_root));
    else
        conds = [];
    end
else
    v_0 = -FF(r(2,2),prime)^(-1)*FF(w_root.^2+r(2,3),prime);
    v_root=v_0.value;
    conds = FF(-r(1,1).*v_root.^2-r(1,2).*v_root-r(1,3)-w_root.^3,prime).value ==0;
end
% Remove roots not satysfying the conditions
v_sol = v_root(conds);
w_sol = w_root(conds);
end

