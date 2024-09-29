function [v_sol,w_sol] = rank3_3C3_235_Fp(r,prime)
%RANK_3C3_235_FP Solve the resulting subsystem of equations for the case R235.
arguments
    r (3,3) {mustBeReal}
    prime (1,1) {mustBeInteger,mustBePositive}
end
% Obtain root candidates from equations
v_root = get_gf_root([-1,-r(1,2),-r(1,3)],prime);
w_root = FF(repmat(-r(3,3),size(v_root)),prime).value;
% Remove roots not satysfying the conditions
conds = FF(-r(2,2).*v_root-r(2,3)-w_root.^2,prime).value ==0;
v_sol = v_root(conds);
w_sol = w_root(conds);
end

