function [params] = get_hyp_cylinder_params(a,b,P,offsets)
%GET_ELL_CYLINDER_PARAMS Summary of this function goes here
%   Detailed explanation goes here
U = -3:0.1:3;
V = -3:0.1:3;
[c1,c2] = meshgrid(U,V);
x =@(s,t,a,b) a.*cosh(s);
y =@(s,t,a,b) b.*sinh(s);
z =@(s,t,a,b) t;
shape = size(c1);
X = x(c1,c2,a,b);
Y = y(c1,c2,a,b);
Z = z(c1,c2,a,b);
tmp = P*[X(:),Y(:),Z(:)]' - P*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
tmp = P*[-X(:),Y(:),Z(:)]' - P*offsets;
XX2 = reshape(tmp(1,:),shape);
YY2 = reshape(tmp(2,:),shape);
ZZ2 = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ},{XX2,YY2,ZZ2}};
end