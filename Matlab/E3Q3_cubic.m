function [result] = E3Q3_cubic(c)
%E3Q3_CUBIC Summary of this function goes here
%   Detailed explanation goes here
arguments
    c (3,20) {mustBeReal} = [1,2,3,-1,4,-2,1,2,0,0,0,0,1,1,1,0,1,1,1,-10;...
        -1,2,3,-1,0,0,1,2,1,0,0,3,0,1,5,-1,7,0,0,-10;...
        2,3,-1,2,0,0,1,1,0,0,1,0,0,1,-3,4,-15,0,0,-10;
        ];
end
clc
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
vars = [x,y,z];
var_vec = [x^3,y^3,z^3,x^2*y,x^2*z,x*y^2,y^2*z,x*z^2,y*z^2,x*y*z,x^2,y^2,z^2,x*y,x*z,y*z,x,y,z,1]';
A = -[c(:,2),c(:,3),c(:,7),c(:,9),c(:,10)*x+c(:,16)]
P = [c(:,6).*x+c(:,12),...
    c(:,8).*x+c(:,13),...
    c(:,4).*x^2+c(:,14).*x+c(:,18),...
    c(:,5).*x^2+c(:,15).*x+c(:,19),...
    c(:,1).*x^3+c(:,11).*x^2+c(:,17).*x+c(:,20)]
Z = pinv(A)*P;
lin_vars = [y^2;z^2;y;z;1]
cube_1 = Z(1,:)*lin_vars;   %y^3
cube_2 = Z(2,:)*lin_vars;   %z^3
quad12 = Z(3,:)*lin_vars;   %y^2*z
quad21 = Z(4,:)*lin_vars;   %y*z^2
mixed = Z(5,:)*lin_vars;    %y*z
% cube_1 = y^3;
% cube_2 = z^3;
% quad12 = y^2*z;
% quad21 = z^2*y;
% mixed = y*z;

p_var = x;
identities = [(cube_1)*lin_vars(4) == (quad12)*lin_vars(3);...      % (y^3*)z = (y^2*z)*y
    (cube_2)*lin_vars(3) == (quad21)*lin_vars(4);...                % (z^3)*y = (z^2*y)*z
    (quad21)*lin_vars(3) == (quad12)*lin_vars(4);...                % (y^2*z)*z = (z^2*y)*y
    (quad12) == (mixed)*lin_vars(3);...                     % (y^2*z) = (y*z)*y
    (mixed)*lin_vars(4) == (quad21)];                      % (y*z)*z = (z^2*y)

sub_new = collect(identities);
old_vars = [lin_vars(3)^3,lin_vars(4)^3,lin_vars(3)^2*lin_vars(4),lin_vars(3)*lin_vars(4)^2];
new_vars = [cube_1,cube_2,quad12,quad21];
sub_new = subs(expand(sub_new),old_vars,new_vars)
expr = lhs(sub_new) - rhs(sub_new);
for i = 1:5
[~,t] = coeffs(expr(i),[x,y,z])
[~,t] = coeffs(expr(i),[y,z])
end
u = sym("u","real");
v = sym("v","real");
w = sym("w","real");

vpa(expr)

sub_vars = [y^2;z^2;y*z;y;z;1]
sub_final = subs(sub_new,sub_vars(1:3),[u;v;w]);
[A_,b] = equationsToMatrix(sub_final,[u;v;w;sub_vars(4:5)]);
A2 = simplify([A_,-b])
size(A2)
A_quad = A2(:,1:3)
A_lin = A2(:,4:6)
pinv(A_quad)
end

