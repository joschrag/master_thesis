function [v_sol,w_sol] = rank3_3C3_125_Fp(r,prime)
%RANK_3C3_125_FP Solve the resulting subsystem of equations for the case R125.
arguments
    r (3,3) {mustBeReal}
    prime (1,1) {mustBeInteger,mustBePositive}
end
% Obtain root candidates from equations
w_root = FF(-r(3,3),prime).value;
v_root = get_gf_root([-1,-r(2,2),-r(2,1)*w_root^2-r(2,3)],prime);
w_root = repmat(w_root,size(v_root));
% Remove roots not satysfying the conditions
conds = FF(-r(1,1).*w_root.^2-r(1,2).*v_root-r(1,3)-w_root.^3,prime).value ==0;
v_sol = v_root(conds);
w_sol = w_root(conds);
end

