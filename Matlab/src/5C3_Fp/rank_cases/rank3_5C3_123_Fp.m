function [v_sol,w_sol] = rank3_5C3_123_Fp(r,prime)
%RANK3_5C3_123_FP Solve the resulting subsystem of equations for the case R123.
arguments
    r (3,2) {mustBeReal}
    prime (1,1) {mustBePrime}
end
% Obtain solution candidates from equations
w_sol = get_gf_root([-1,-r(2,:)],prime);
if isempty(w_sol)
    v_sol = [];
    return
end
v_sol = FF(-(r(3,1).*w_sol+r(3,2)),prime).value;
% Remove roots not satysfying the conditions
control = FF(-(r(1,1).*w_sol+v_sol.^2+r(1,2)),prime).value == 0;
w_sol = w_sol(control);
v_sol = v_sol(control);
end