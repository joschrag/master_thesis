function [result] = E3Q3(c)
%E3Q3 Summary of this function goes here
%   Detailed explanation goes here
arguments
c (3,10) {mustBeReal} = [1,2,1,1,1,0,1,1,1,-10;...
    1,1,3,0,1,-1,7,0,0,-10;...
    1,1,1,0,0,1,-15,0,0,-10;
    ];
end
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");

A = -[c(:,2),c(:,3),c(:,6)];
rref(A)
P = [c(:,4).*x+c(:,8),c(:,5).*x+c(:,9),c(:,1).*x^2+c(:,7).*x+c(:,10)];
if rank(A) ~= 3
    P2 = pinv(A)*P; %idk
else
    P2 = A\P;
end

lin_vars = [y;z;1];
y_quad = P2(1,:)*lin_vars;
z_quad = P2(2,:)*lin_vars;
yz = P2(3,:)*lin_vars;

identities = [(y_quad)*z == (yz)*y;...
    (yz)*z == (z_quad)*y;...
    (yz)*(yz) == y_quad*z_quad];

sub1 = collect(identities);
sub2 = subs(sub1,[y^2,z^2,y*z],[y_quad,z_quad,yz]);
sub3 = subs(expand(sub2),[y^2,z^2,y*z],[y_quad,z_quad,yz]);
[A_,b] = equationsToMatrix(sub3,[y,z]);
A2 = [A_,-b];
p = det(A2);
coef = coeffs(p,"All");
p_root = sturm(eval(coef));
if size(p_root,1) == 0
    fprintf("No common intersection points found\n")
    plot_and_color(c)
    return
end
extra = abs(max(p_root)-min(p_root))*0.35;
ex_x_roots = [min(p_root)-extra,max(p_root)+extra];
result = zeros(numel(p_root),3);
result(:,1) = p_root;
for i =1:numel(p_root)
    M = subs(A2,x,p_root(i));
    [~,S,V] = svd(M);
    [~,index] = min(diag(S));
    result(i,2) = V(1,index)/V(3,index);
    result(i,3) = V(2,index)/V(3,index);    
end
ax = gca();
scatter3(ax,result(:,1),result(:,2),result(:,3),25,"k","filled")
hold(ax,"on")
min_y = min(result(:,2));
max_y = max(result(:,2));
min_z = min(result(:,3));
max_z = max(result(:,3));
extra = abs(max_y-min_y)*0.35;
ex_y_roots = [min_y-extra,max_y+extra];
extra = abs(max_z-min_z)*0.35;
ex_z_roots = [min_z-extra,max_z+extra];
plot_and_color(c)
axis(horzcat(ex_x_roots,ex_y_roots,ex_z_roots))
end

