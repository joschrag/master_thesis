function [v_sol,w_sol] = rank3_5C3_124_Fp(r,prime)
%RANK3_5C3_124_FP Solve the resulting subsystem of equations for the case R124.
arguments
    r (3,2) {mustBeReal}
    prime (1,1) {mustBePositive,mustBeInteger}
end
% Obtain solutions from equations
v_root = get_gf_root([-1,-r(1,1),-r(1,2)],prime);
w_root = FF(repmat(-r(3,2),size(v_root)),prime).value;
% Remove roots not satysfying the conditions
control = FF(-w_root.^2-r(2,1).*v_root-r(2,2),prime).value == 0;
v_sol = v_root(control);
w_sol = w_root(control);
end