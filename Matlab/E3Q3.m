function [] = E3Q3()
%E3Q3 Summary of this function goes here
%   Detailed explanation goes here

x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
vars = [x^2,y^2,z^2,x*y,x*z,y*z,x,y,z,1];
c = [4,0,0,1,0,1,1,1,1,-1;...
          -3,1,1,0,1,0,7,0,0,-1;...
          2,-1,2,1,0,1,-15,0,0,-1;
         ];
A = -[c(:,2),c(:,3),c(:,6)];
P = [c(:,4).*x+c(:,8),c(:,5).*x+c(:,9),c(:,1).*x^2+c(:,7).*x+c(:,10)];
rref(inv(A))
P2 = A\P;
lin_vars = [y;z;1];
y_quad = P2(1,:)*lin_vars;
z_quad = P2(2,:)*lin_vars;
yz = P2(3,:)*lin_vars;


identities = [(y_quad)*z == (yz)*y;...
    (yz)*z == (z_quad)*y;...
    (yz)*(yz) == y_quad*z_quad];

sub1 = simplify(expand(identities));
sub2 = subs(sub1,[y^2,z^2,y*z],[y_quad,z_quad,yz])
[A_,b] = equationsToMatrix(sub2,[y,z])
A2 = simplify([A_,-b])
det(A2)



end

