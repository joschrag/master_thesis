function [u_sol,v_sol] = rank2_5C3_12(r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (2,3) {mustBeReal}
end
u = sym("u","real");
v = sym("v","real");
u_sol = [];
v_sol = [];
if r(1,2) ~= 0
    v_0 = -(u^2+r(1,1)*u+r(1,3))/r(1,2);
    pol_2 = -subs(r(2,1)*u+r(2,2)*v+v^2+r(2,3),v,v_0);
    u_sol = roots(coeffs(pol_2,"All"));
    u_sol = u_sol(imag(u_sol)==0);
    v_sol = subs(v_0,u,u_sol);
    
elseif r(2,1) ~= 0
    t_0 = -(v^2+r(2,2)*v+r(2,3))/r(2,1);
    pol_1 = -subs(r(1,1)*u+r(1,2)*v+u^2+r(1,3),u,t_0);
    v_sol_2 = roots(coeffs(pol_1,"All"));
    v_sol_2 = v_sol_2(imag(v_sol_2)==0);
    v_sol = [v_sol;v_sol_2];
    u_sol = [u_sol;subs(t_0,v,v_sol_2)];
end
if r(1,2) == 0 && r(2,1) == 0
    t_0 = roots([-1,-r(1,1),-r(1,3)]);
    t_0 = t_0(imag(t_0)==0);
    v_0 = roots([-1,-r(2,2),-r(2,3)]);
    v_0 = v_0(imag(v_0)==0);
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