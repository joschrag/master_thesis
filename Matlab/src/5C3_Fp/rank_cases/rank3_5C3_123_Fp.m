function [u_sol,v_sol] = rank3_5C3_123_Fp(r,prime)
%RANK3_5C3_123_FP Solve the resulting subsystem of equations for the case R123.
arguments
    r (3,2) {mustBeReal}
    prime (1,1) {mustBePrime}
end
% Obtain solution candidates from equations
v_sol = get_gf_root([-1,-r(2,:)],prime);
if isempty(v_sol)
    u_sol = [];
    return
end
u_sol = FF(-(r(3,1).*v_sol+r(3,2)),prime).value;
% Remove roots not satysfying the conditions
control = FF(-(r(1,1).*v_sol+u_sol.^2+r(1,2)),prime).value == 0;
v_sol = v_sol(control);
u_sol = u_sol(control);
end