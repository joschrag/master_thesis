function [v_sol,w_sol] = rank2_5C3_12(r)
%RANK2_5C3_12 Solve the resulting subsystem of equations for the case R12.
arguments
    r (2,3) {mustBeReal}
end
v = sym("v","real");
w = sym("w","real");
v_sol = [];
w_sol = [];
% Obtain solutions from equation based on coefficients
if r(1,2) ~= 0
    w_0 = -(v^2+r(1,1)*v+r(1,3))/r(1,2);
    pol_2 = -subs(r(2,1)*v+r(2,2)*w+w^2+r(2,3),w,w_0);
    v_sol = roots(coeffs(pol_2,"All"));
    v_sol = v_sol(imag(v_sol)==0);
    w_sol = subs(w_0,v,v_sol);
elseif r(2,1) ~= 0
    v_0 = -(w^2+r(2,2)*w+r(2,3))/r(2,1);
    pol_1 = -subs(r(1,1)*v+r(1,2)*w+v^2+r(1,3),v,v_0);
    w_sol_2 = roots(coeffs(pol_1,"All"));
    w_sol_2 = w_sol_2(imag(w_sol_2)==0);
    w_sol = [w_sol;w_sol_2];
    v_sol = [v_sol;subs(v_0,w,w_sol_2)];
end
if r(1,2) == 0 && r(2,1) == 0
    v_0 = roots([-1,-r(1,1),-r(1,3)]);
    v_0 = v_0(imag(v_0)==0);
    w_0 = roots([-1,-r(2,2),-r(2,3)]);
    w_0 = w_0(imag(w_0)==0);
    if isempty(v_0) || isempty(w_0)
        w_sol = [];
        v_sol = [];
        return
    end
    [v_0,w_0] = meshgrid(v_0,w_0);
    v_sol = [v_sol;v_0];
    w_sol = [w_sol;w_0];
end
single_sols = unique([w_sol,v_sol],"rows");
if isempty(single_sols)
    w_sol = [];
    v_sol = [];
else
    w_sol = [single_sols(:,1)];
    v_sol = [single_sols(:,2)];
end
end