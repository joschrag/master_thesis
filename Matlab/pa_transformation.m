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
ev_count = {pos_eig_vals,neg_eig_vals,zero_eig_vals};
res = classify_quadrics(ev_count,const,g,h);

switch res
    case 3
        params = get_ellipsoid_params(a,b,c,V,offsets);
        dname = "ellipsoid";
    case 4
        params = get_hyperboloid_1_params(a,b,c,V,offsets);
        dname = "hyperboloid1";
    case 5
        params = get_hyperboloid_2_params(a,b,c,V,offsets);
        dname = "hyperboloid2";
    case 6
        params = get_ell_paraboloid_params(a,b,V,offsets);
        dname = "ell paraboloid";
    case 7
        params = get_hyp_paraboloid_params(a,b,V,offsets);
        dname = "hyp paraboloid";
    case 2
        params = get_ell_cone_params(a,b,c,V,offsets);
        dname = "ell cone";
    case 8
        params = get_ell_cylinder_params(a,b,V,offsets);
        dname = "ell cylinder";
    case 10
        params = get_hyp_cylinder_params(a,b,V,offsets);
        dname = "hyp cylinder";
    case 11
        tmp = [g,h];
        a = eval(sqrt(abs(tmp(tmp~=0)/(2*eig_vals(eig_vals~=0)))))
        params = get_par_cylinder_params(a,V,offsets);
        dname = "par cylinder";
    case 12
%         fun2(V*new_vars) + fun1(V*new_vars) + const
        a = sqrt(eig_vals(eig_vals~=0)./(-const));
        params = get_parallel_plane_params(a, V, offsets);
        dname = "par planes";
    case 13
%         fun2(V*new_vars) + fun1(V*new_vars) + const
        params = get_plane_params(a, b, V, offsets);
        dname = "plane";
    case 14
%         fun2(V*new_vars) + fun1(V*new_vars) + const
        tmp = sqrt(abs(eig_vals(eig_vals~=0)));
        a = tmp(1);
        b = tmp(2);
        params = get_cross_plane_params(a,b,V,offsets);
        dname = "cross planes";
    case 9
%         fun2(V*new_vars) + fun1(V*new_vars) + const
        params = get_line_params(V,offsets);
        dname = "line";
    case 1
        params = {{0,0,0}};
        dname = "point";
end
ax = gca();
for i = 1:numel(params)
    coords = params{i};
    if all(size(coords{3})>1)
    go = surf(ax,coords{1},coords{2},coords{3});
    set(go,"EdgeColor","None")
    elseif any(size(coords{3})>1)
        go = plot3(ax,coords{1},coords{2},coords{3});
    else
        go = scatter3(ax,coords{1},coords{2},coords{3});
    end
    go.DisplayName = dname;
    hold(ax,"on")
    disp(fun2(V*new_vars) + fun1(V*new_vars) + const)
end


function [g, h] = detect_linear_factors(new_vars, fun1, V, fun2)
g = 0;
h = 0;
single_vars = new_vars(~ismember(new_vars,symvar(fun1(V*new_vars))));
if numel(single_vars) >= 1
    tmp = coeffs(fun2(V*new_vars),single_vars(1),"All");
    if numel(tmp) == 2
        g = tmp(1);
    end
end
if numel(single_vars) >= 2
    tmp = coeffs(fun2(V*new_vars),single_vars(2),"All");
    if numel(tmp) == 2
        h = tmp(1);
    end
end