function [result] = E3Q3(c)
%E3Q3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    c (:,10) {mustBeReal} = [...
        1,2,1,1,1,0,1,1,1,-10;...
        1,1,3,0,1,-1,7,0,0,-10;...
        1,1,1,0,0,1,-15,0,0,-10;...
        ];
end
clc;
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
vars = [x,y,z];
v = [x^2,y^2,z^2,x*y,x*z,y*z,x,y,z,1]';
A = -[c(:,2),c(:,3),c(:,6)];
B = -[c(:,1),c(:,3),c(:,5)];
C = -[c(:,1),c(:,2),c(:,4)];
[m,~] = size(c);
r_A = rank(A);
r_B = rank(B);
r_C = rank(C);
if any([r_A,r_B,r_C] == 3)
    if m == 3
        if r_A == 3
            fprintf("Using P(x)\n")
            Q = A;
            P = [c(:,4).*x+c(:,8),c(:,5).*x+c(:,9),c(:,1).*x^2+c(:,7).*x+c(:,10)];
            lin_vars = [y;z;1];
            p_var = x;
        elseif r_B == 3
            fprintf("Using P(y)\n")
            Q = B;
            P = [c(:,4).*y+c(:,7),c(:,6).*y+c(:,9),c(:,2).*y^2+c(:,8).*y+c(:,10)];
            lin_vars = [x;z;1];
            p_var = y;
        elseif r_C == 3
            fprintf("Using P(z)\n")
            Q = C;
            P = [c(:,5).*z+c(:,7),c(:,6).*z+c(:,8),c(:,3).*z^2+c(:,9).*z+c(:,10)];
            lin_vars = [x;y;1];
            p_var = z;
        end
        P2 = Q\P;
    else
        fprintf("Using MP-inverse\n")
        P_X = [c(:,4).*x+c(:,8),c(:,5).*x+c(:,9),c(:,1).*x^2+c(:,7).*x+c(:,10)];
        P2 = pinv(A)*P_X; %idk
        lin_vars = [y;z;1];
        p_var = x;
        pinv(A)*A
    end
else
    fprintf("Using MP-inverse\n")
    P_X = [c(:,4).*x+c(:,8),c(:,5).*x+c(:,9),c(:,1).*x^2+c(:,7).*x+c(:,10)];
    P2 = pinv(A)*P_X; %idk
    lin_vars = [y;z;1];
    p_var = x;
    pinv(A)*A
end

quad_1 = P2(1,:)*lin_vars
quad_2 = P2(2,:)*lin_vars
mixed = P2(3,:)*lin_vars

identities = [(quad_1)*lin_vars(2) == (mixed)*lin_vars(1);...
    (mixed)*lin_vars(2) == (quad_2)*lin_vars(1);...
    (mixed)*(mixed) == quad_1*quad_2];

sub1 = collect(identities);
old_vars = [lin_vars(1)^2,lin_vars(2)^2,lin_vars(1)*lin_vars(2)];
new_vars = [quad_1,quad_2,mixed];
sub2 = subs(sub1,old_vars,new_vars);
sub3 = subs(expand(sub2),old_vars,new_vars);
vpa(collect(simplify(lhs(sub3)-rhs(sub3)),[x,y,z] ),4)
[A_,b] = equationsToMatrix(sub3,lin_vars(1:2));
A2 = [A_,-b];
vpa(collect(A2))
p = det(A2);
coef = coeffs(p,p_var,"All");
vpa(coef)
p_root = sturm(eval(coef))
idx = [find(vars==p_var),find(vars==lin_vars(1)),find(vars==lin_vars(2))];
if size(p_root,1) == 0
    plot_and_color(c)
    error("No common intersection points found\n")
end
result = zeros(numel(p_root),3);
result(:,idx(1)) = p_root;
for i =1:numel(p_root)
    M = subs(A2,p_var,p_root(i))
    [~,S,V] = svd(M)
    [~,index] = min(diag(S));
    result(i,idx(2)) = V(1,index)/V(3,index);
    result(i,idx(3)) = V(2,index)/V(3,index);
end
ax = gca();
s = scatter3(ax,result(:,1),result(:,2),result(:,3),25,"k","filled");
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
n = zeros(numel(p_root),1);
for i=1:size(result,1)
    n(i) = norm(eval(subs(c*v,[x,y,z],result(i,:))));
    eval_str = join(repmat("\n%g",[1,m]),"");
    line_break = join(repmat("=",[1,10]),"");
    result_str = sprintf("Lösung %%i -> \nNorm: %%g%s\n%%s\n",eval_str);
    fprintf(result_str,i,n(i),eval(subs(c*v,[x,y,z],result(i,:))),line_break);
end
end

