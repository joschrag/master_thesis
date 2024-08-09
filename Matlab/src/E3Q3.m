function [result] = E3Q3(c)
%E3Q3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    c (:,10) {mustBeReal} = c = [...
        1,1,1,2,0,1,1,1,1,-10;...
        1,0,1,1,-1,3,7,0,0,-10;...
        1,0,0,1,1,1,-15,0,0,-10;...
        ];
end
clc;
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
vars = [x,y,z];
v = [x^2,x*y,x*z,y^2,y*z,z^2,x,y,z,1]';
[m,~] = size(c);
[lin_vars, p_var, P2] = split_matrices(c,m, x, y, z);

quad_1 = P2(1,:)*lin_vars;
quad_2 = P2(3,:)*lin_vars;
mixed = P2(2,:)*lin_vars;

identities = [(quad_1)*lin_vars(2) == (mixed)*lin_vars(1);...
    (mixed)*lin_vars(2) == (quad_2)*lin_vars(1);...
    (mixed)*(mixed) == quad_1*quad_2];

sub1 = collect(identities);
old_vars = [lin_vars(1)^2,lin_vars(2)^2,lin_vars(1)*lin_vars(2)];
new_vars = [quad_1,quad_2,mixed];
sub2 = subs(sub1,old_vars,new_vars);
sub3 = subs(expand(sub2),old_vars,new_vars);
[A_,b] = equationsToMatrix(sub3,lin_vars(1:2));
A2 = [A_,-b];
p = det(A2);
coef = coeffs(p,p_var,"All");
p_root = sturm(double(coef));
% p_root = roots(coef);
% p_root = p_root(abs(imag(p_root))<10^(-10));
idx = [find(vars==p_var),find(vars==lin_vars(1)),find(vars==lin_vars(2))];
if size(p_root,1) == 0
    plot_and_color(c)
    error("No common intersection points found")
end
result = zeros(numel(p_root),3);
result(:,idx(1)) = p_root;
for i =1:numel(p_root)
    M = subs(A2,p_var,p_root(i));
    %vpa(subs(det(A2),x,p_root(i)))
    [~,S,V] = svd(M);
    [~,index] = min(diag(S));
    result(i,idx(2)) = V(1,index)/V(3,index);
    result(i,idx(3)) = V(2,index)/V(3,index);
end
ax = gca();
s = scatter3(ax,result(:,1),result(:,2),result(:,3),30,"k","filled");
set(s,"DisplayName","solution")
hold(ax,"on")
min_x = min(result(:,1));
max_x = max(result(:,1));
min_y = min(result(:,2));
max_y = max(result(:,2));
min_z = min(result(:,3));
max_z = max(result(:,3));
min_dist = 0.1;
extra = abs(max(max_x-min_x,min_dist))*0.35;
ex_x_roots = [min_x-extra,max_x+extra];
extra = abs(max(max_y-min_y,min_dist))*0.35;
ex_y_roots = [min_y-extra,max_y+extra];
extra = abs(max(max_z-min_z,min_dist))*0.35;
ex_z_roots = [min_z-extra,max_z+extra];
plot_and_color(c)
axis(horzcat(ex_x_roots,ex_y_roots,ex_z_roots))
equations = c*v;
if numel(result) > 0
    print_solutions(result,equations,x,y,z)
end
end

