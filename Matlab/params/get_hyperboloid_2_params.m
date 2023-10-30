function [params] = get_hyperboloid_2_params(a,b,c,V,offsets)
%GET_HYPERBOLOID_2_PARAMS Summary of this function goes here
%   Detailed explanation goes here
step_arc = 0.01;
x =@(s,t,a,b,c) a.*sinh(s).*cos(t);
y =@(s,t,a,b,c) b.*sinh(s).*sin(t);
z =@(s,t,a,b,c) c.*cosh(s);
[c1,c2] = meshgrid(-2:0.1:2,0:step_arc:2*pi+step_arc);
shape = size(c1);
X = x(c1,c2,a,b,c);
Y = y(c1,c2,a,b,c);
Z = z(c1,c2,a,b,c);
tmp = V*[X(:),Y(:),Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
tmp = V*[X(:),Y(:),-Z(:)]' - V*offsets;
XX2 = reshape(tmp(1,:),shape);
YY2 = reshape(tmp(2,:),shape);
ZZ2 = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ},{XX2,YY2,ZZ2}};
end