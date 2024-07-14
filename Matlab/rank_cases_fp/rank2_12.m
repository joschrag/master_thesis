function [u_sol,v_sol] = rank2_12(r,prime)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (2,3) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePositive,mustBeInteger}
end
u = sym("t",["real","integer"]);
v = sym("u",["real","integer"]);
u_sol = [];
v_sol = [];
if r(1,2) ~= 0
    v_0 = FF(-(u^2+r(1,1)*u+r(1,3)),prime,u)*FF(r(1,2),prime)^(-1);
    pol_2 = FF(-subs(r(2,1)*u+r(2,2)*v+v^2+r(2,3),v,v_0.value),prime);
    u_sol = get_gf_root(coeffs(pol_2.value,"All"),prime);
    v_sol = subs(v_0,u,u_sol).value;
elseif r(2,1) ~= 0
    t_0 = FF(-(v^2+r(2,2)*v+r(2,3)),prime)*FF(r(2,1),prime)^(-1);
    pol_1 = FF(-subs(r(1,1)*u+r(1,2)*v+u^2+r(1,3),u,t_0.value),prime);
    v_sol_2 = get_gf_root(coeffs(pol_1.value,"All"),prime);
    v_sol = [v_sol;v_sol_2];
    u_sol = [u_sol;subs(t_0,v,v_sol_2).value];
end
if r(1,2) == 0 && r(2,1) == 0
    t_0 = get_gf_root([-1,-r(1,1),-r(1,3)],prime);
    v_0 = get_gf_root([-1,-r(2,2),-r(2,3)],prime);
    if isempty(t_0) || isempty(v_0)
        v_sol = [];
        u_sol = [];
        return
    end
    tmp = zeros(numel(v_0)*numel(t_0),2);
    for i=1:numel(t_0)
        for j=1:numel(v_0)
            tmp((i-1)*numel(v_0)+j,:) = [t_0(i),v_0(j)];
        end
    end
    u_sol = [u_sol;tmp(:,1)];
    v_sol = [v_sol;tmp(:,2)];
end

single_sols = unique([v_sol,u_sol],"rows");
if isempty(single_sols)
    v_sol = [];
    u_sol = [];
else
    v_sol = [single_sols(:,1)];
    u_sol = [single_sols(:,2)];
end
end