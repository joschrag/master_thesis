function [v_sol,w_sol] = rank4_3C3_1234_Fp(r,prime)
%RANK_3C3_1234_FP Solve the resulting subsystem of equations for the case R1234.
arguments
    r (4,2) {mustBeReal}
    prime (1,1) {mustBePrime}
end
% Obtain root candidates from equations
w_root = get_gf_root([-1,-r(3,1),-r(3,2)],prime);
v_root = FF(-r(4,1).*w_root-r(4,2),prime).value;
% Remove roots not satysfying the conditions
conds = FF(-r(1,1).*w_root-r(1,2)-w_root.^3,prime).value == 0 & FF(-r(2,1).*w_root-r(2,2)-v_root.^2,prime).value == 0;
v_sol = v_root(conds);
w_sol = w_root(conds);
end

