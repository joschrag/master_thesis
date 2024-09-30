function [u_sol,v_sol] = rank3_5C3_124_Fp(r,prime)
%RANK3_5C3_124_FP Solve the resulting subsystem of equations for the case R124.
arguments
    r (3,2) {mustBeReal}
    prime (1,1) {mustBePositive,mustBeInteger}
end
% Obtain solutions from equations
u_0 = get_gf_root([-1,-r(1,1),-r(1,2)],prime);
v_0 = FF(repmat(-r(3,2),size(u0)),prime).value;
% Remove roots not satysfying the conditions
control = FF(-v_0^2-r(2,1).*u_0-r(2,2),prime).value == 0;
u_sol = u_0(control);
v_sol = v_0(control);
end