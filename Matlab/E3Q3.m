function [] = E3Q3()
%E3Q3 Summary of this function goes here
%   Detailed explanation goes here

x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
variables = [x^2;y^2;z^2;x*y;x*z;y*z;x;y;z;1];
c = [1,2,1,1,0,1,1,1,1,-10;...
    1,6,1,0,1,-1,7,0,0,-10;...
    1,-1,-1,1,0,1,-15,0,0,-10;
    ];
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

sub1 = simplify(collect(identities));
sub2 = subs(sub1,[y^2,z^2,y*z],[y_quad,z_quad,yz]);
sub3 = subs(expand(sub2),[y^2,z^2,y*z],[y_quad,z_quad,yz]);
[A_,b] = equationsToMatrix(sub3,[y,z]);
A2 = simplify([A_,-b]);
p = det(A2);
coef = coeffs(p,"All");
p_root = sturm(double(coef));
extra = abs(max(p_root)-min(p_root))*0.35;
ex_x_roots = [min(p_root)-extra,max(p_root)+extra];
ex_y_roots = [inf,-inf];
ex_z_roots = [inf,-inf];
for x_root = p_root'
    M = subs(A2,x,x_root);
    %eigenvector approach
    [~,S,V] = svd(M);
    [~,index] = min(diag(S));
    y_root = V(1,index)/V(3,index);
    z_root = V(2,index)/V(3,index);
    ex_y_roots = [min(ex_y_roots(1),y_root),max(ex_y_roots(2),y_root)];
    ex_z_roots = [min(ex_z_roots(1),z_root),max(ex_z_roots(2),z_root)];
    ax = gca();
    scatter3(ax,x_root,y_root,z_root,25,"k","filled")
    hold(ax,"on")
    double(subs(c*variables,[x,y,z],[x_root,y_root,z_root]))
end
extra = abs(ex_y_roots(2)-ex_y_roots(1))*0.35;
ex_y_roots = [ex_y_roots(1)-extra,ex_y_roots(2)+extra];
extra = abs(ex_z_roots(2)-ex_z_roots(1))*0.35;
ex_z_roots = [ex_z_roots(1)-extra,ex_z_roots(2)+extra];
plot_and_color(c)
axis(eval(horzcat(ex_x_roots,ex_y_roots,ex_z_roots)))
end

