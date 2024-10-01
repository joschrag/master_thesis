function [v_sol,w_sol] = rank4_4C3_1235_Fp(r,prime)
%RANK4_4C3_1235_FP Solve the resulting subsystem of equations for the case R1235.
arguments
    r (4,2) {mustBeReal}
    prime (1,1) {mustBePrime}
end
% Obtain solution candidates from equations
v_root = get_gf_root([-1,-r(1,1),-r(1,2)],prime);
w_root = FF(repmat(-r(4,2),size(v_root)),prime).value;
% Remove roots not satysfying the conditions
conds = FF(-r(3,1).*v_root-r(3,2)-w_root.^2,prime).value == 0 & FF(-r(2,1).*v_root-r(2,2)-v_root.*w_root,prime).value == 0;
v_sol = v_root(conds);
w_sol = w_root(conds);
end