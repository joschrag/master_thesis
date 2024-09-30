function [v_sol,w_sol] = rank2_5C3_23_Fp(r,prime)
%RANK2_5C3_23_FP Solve the resulting subsystem of equations for the case R23.
arguments
    r (2,3) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePrime}
end
% Obtain solutions from equations
w_0 = get_gf_root([-1,-r(1,2),-r(1,3)],prime);
if ~isempty(w_0)
    w_sol = w_0;
    v_sol = FF(-r(2,2).*w_0-r(2,3),prime).value;
else
    w_sol = [];
    v_sol = [];
end
end