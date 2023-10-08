function pa_transformation(coeff_vec)
%PA_TRANSFORMATION Summary of this function goes here
%   Detailed explanation goes here
%   arguments:
%       coeffs: coefficient vector of (x^2,y^2,z^2,xy,xz,yz,x,y,z,1)
arguments
    coeff_vec (1,10) {mustBeReal} = ones(1,10);
end
% if abs(coeff_vec(10)) >= tol
%     coeff_vec = coeff_vec./(-coeff_vec(10));
%     rhs = 1;
% else
%     rhs = 0;
% end
A = [coeff_vec(1),coeff_vec(4)/2,coeff_vec(5)/2;...
    coeff_vec(4)/2,coeff_vec(2),coeff_vec(6)/2;...
    coeff_vec(5)/2,coeff_vec(6)/2,coeff_vec(3)];
B = [coeff_vec(7),coeff_vec(8),coeff_vec(9)];
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
var = [x,y,z,1];
fun1 = matlabFunction(simplify(expand(var(1:3)*A*var(1:3)')),'Vars',{[x;y;z]});
fun2 = matlabFunction(simplify(expand(dot(B,var(1:3)))),'Vars',{[x;y;z]});

[V,D] = eig(A);
V = [V(:,1)./norm(V(:,1)),V(:,2)./norm(V(:,2)),V(:,3)./norm(V(:,3))];
eig_vals = diag(D);

if round(det(V)) == -1
    V = [V(:,2),V(:,1),V(:,3)];
    eig_vals = [eig_vals(2),eig_vals(1),eig_vals(3)];
end
lambda1 = eig_vals(1);
lambda2 = eig_vals(2);
lambda3 = eig_vals(3);
pos_eig_vals = sum(eig_vals > 0);
neg_eig_vals = sum(eig_vals < 0);
zero_eig_vals = 3-pos_eig_vals-neg_eig_vals;
xi = sym("xi","real");
eta = sym("eta","real");
zeta = sym("zeta","real");
new_vars = [xi;eta;zeta];
const = coeff_vec(10);
offsets = zeros(3,1);
if any(B~=0)
    lhs_lin = simplify(fun2(V*new_vars));
else
    lhs_lin = 0;
end
if lambda1 ~= 0
    [tmp1,~] = coeffs(lhs_lin,xi,"All");
    if numel(tmp1) > 1
        factor = double(tmp1(1));
        [offsets(1),const] = square_complete(lambda1,factor,const);
    else
        offsets(1) = 0;
    end
end
if lambda2 ~= 0
    [tmp1,~] = coeffs(lhs_lin,eta,"All");
    if numel(tmp1) > 1
        factor = double(tmp1(1));
        [offsets(2),const] = square_complete(lambda2,factor,const);
    else
        offsets(2) = 0;
    end
end
if lambda3 ~= 0
    [tmp1,~] = coeffs(lhs_lin,zeta,"All");
    if numel(tmp1) > 1
        factor = double(tmp1(1));
        [offsets(3),const] = square_complete(lambda3,factor,const);
    else
        offsets(3) = 0;
    end

end
if const ~= 0
    tmp = 1./(sqrt(abs(eig_vals./(-const))));
else
    tmp = 1./(sqrt(abs(eig_vals)));
end
if lambda1 ~= 0
    a = double(tmp(1));
else
    a = 0;
end
if lambda2 ~= 0
    b = double(tmp(2));
else
    b = 0;
end
if lambda3 ~= 0
    c = double(tmp(3));
else
    c = 0;
end
[g, h] = detect_linear_factors(new_vars, fun1, V, fun2);
ev_count = {pos_eig_vals,neg_eig_vals,zero_eig_vals}
res = classify_quadrics(ev_count,const,g,h)

switch res
    case "ellipsoid"
        plot_ellipsoid(a,b,c,V,offsets)
    case "hyperboloid 1"
        plot_hyperboloid_1(a,b,c,V,offsets)
    case "hyperboloid 2"
        plot_hyperboloid_2(a,b,c,V,offsets)
    case "ell paraboloid"
        plot_ell_paraboloid(a,b,c,V,offsets)
    case "hyp paraboloid"
        plot_hyp_paraboloid(a,b,c,V,offsets)
    case "ell cone"
        plot_ell_cone(a,b,c,V,offsets)
    case "ell cylinder"
        plot_ell_cylinder(a,b,c,V,offsets)
    case "hyp cylinder"
        plot_hyp_cylinder(a,b,c,V,offsets)
    case "par cylinder"
        plot_par_cylinder(a,b,c,V,offsets)
    case "two planes"
        fun2(V*new_vars) + fun1(V*new_vars) + const
    case "one plane"
        fun2(V*new_vars) + fun1(V*new_vars) + const
    case "one line"
        fun2(V*new_vars) + fun1(V*new_vars) + const
    case "point solution"
        fun2(V*new_vars) + fun1(V*new_vars) + const
end








if zero_eig_vals == 0
    if const ~= 0 % cases [+++1,++-1,+--1]
        plot_quadric_3d(pos_eig_vals,neg_eig_vals,a,b,c,V,offsets);
    else % cases [++-0]
        plot_ell_cone(a,b,c,V,offsets)
    end
elseif zero_eig_vals == 1
    if const == 0
        if pos_eig_vals == 1 % case [+-00]
            disp("hyp Para")
            single_var = new_vars(~ismember(new_vars,symvar(fun1(V*new_vars))));
            cc = coeffs(fun2(V*new_vars),single_var,"All");
            tmp_fac = double(2.*eig_vals(eig_vals~=0)./(sign(cc(1))*cc(1)));
            a = abs(tmp_fac(1));
            b = abs(tmp_fac(2));
            plot_paraboloids(pos_eig_vals,a,b,V,offsets)
        elseif pos_eig_vals == 2 % case [++00]
            disp("ell Para")
            single_var = new_vars(~ismember(new_vars,symvar(fun1(V*new_vars))));
            cc = coeffs(fun2(V*new_vars),single_var,"All");
            tmp_fac = double(2.*eig_vals(eig_vals~=0)./(sign(cc(1))*cc(1)));
            a = tmp_fac(1);
            b = tmp_fac(2);
            plot_paraboloids(pos_eig_vals,a,b,V,offsets)
        end
    else % case [+-01,++01]
        tmp = [a,b,c];
        tmp = tmp(tmp~=0);
        a = tmp(1);
        b = tmp(2);
        plot_cylinder(pos_eig_vals,a,b,V,offsets)

        if pos_eig_vals == 1
            disp("ell Zyl")

        elseif pos_eig_vals == 2
            disp("hyp Zyl")
        end
    end
elseif zero_eig_vals == 2 % cases [+000,-000]
    single_var = new_vars(ismember(new_vars,symvar(fun2(V*new_vars))));
    cc = coeffs(fun2(V*new_vars),single_var,"All");
    tmp_fac = double(2.*eig_vals(eig_vals~=0)./(sign(cc(1))*cc(1)));
    a = tmp_fac(1);
    plot_parabola(a,V,offsets)
end

function [g, h] = detect_linear_factors(new_vars, fun1, V, fun2)
g = 0;
h = 0;
single_vars = new_vars(~ismember(new_vars,symvar(fun1(V*new_vars))));
tmp = coeffs(fun2(V*new_vars),single_vars(1),"All");
if numel(tmp) == 2
    g= tmp(1);
end
if numel(single_vars) == 2
    tmp = coeffs(fun2(V*new_vars),single_vars(2),"All");
    if numel(tmp) == 2
        h = tmp(1);
    end
end