function [v_sol,w_sol] = rank2_5C3_12_Fp(r,prime)
%RANK2_5C3_12_FP Solve the resulting subsystem of equations for the case R12.
arguments
    r (2,3) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePrime}
end
v = sym("v",["real","integer"]);
w = sym("w",["real","integer"]);
v_sol = [];
w_sol = [];
% Obtain results from equations
if r(1,2) ~= 0
    w_0 = FF(-(v^2+r(1,1)*v+r(1,3)),prime,v)*FF(r(1,2),prime)^(-1);
    pol_2 = FF(-subs(r(2,1)*v+r(2,2)*w+w^2+r(2,3),w,w_0.value),prime);
    v_sol = get_gf_root(coeffs(pol_2.value,"All"),prime);
    w_sol = subs(w_0,v,v_sol).value;
elseif r(2,1) ~= 0
    v_0 = FF(-(w^2+r(2,2)*w+r(2,3)),prime)*FF(r(2,1),prime)^(-1);
    pol_1 = FF(-subs(r(1,1)*v+r(1,2)*w+v^2+r(1,3),v,v_0.value),prime);
    w_sol_2 = get_gf_root(coeffs(pol_1.value,"All"),prime);
    w_sol = [w_sol;w_sol_2];
    v_sol = [v_sol;subs(v_0,w,w_sol_2).value];
end
if r(1,2) == 0 && r(2,1) == 0
    v_0 = get_gf_root([-1,-r(1,1),-r(1,3)],prime);
    w_0 = get_gf_root([-1,-r(2,2),-r(2,3)],prime);
    [v_0,w_0] = meshgrid(v_0,w_0);
    v_sol = [v_sol;v_0];
    w_sol = [w_sol;w_0];
end
% Remove duplicate solutions
single_sols = unique([w_sol,v_sol],"rows");
if isempty(single_sols)
    w_sol = [];
    v_sol = [];
else
    w_sol = [single_sols(:,1)];
    v_sol = [single_sols(:,2)];
end
end