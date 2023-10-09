function [params] = get_hyperboloid_1_params(a,b,c,V,offsets)
%GET_HYPERBOLOID_1_PARAMS Summary of this function goes here
%   Detailed explanation goes here
step_arc = 0.01;
x =@(s,t,a,b,c)  a.*cosh(s).*cos(t);
y =@(s,t,a,b,c) b.*cosh(s).*sin(t);
z =@(s,t,a,b,c) c.*sinh(s);
[c1,c2] = meshgrid(0:step_arc:2*pi+step_arc,-10:0.1:10);
shape = size(c1);
X = x(c1,c2,a,b,c);
Y = y(c1,c2,a,b,c);
Z = z(c1,c2,a,b,c);
tmp = V*[X(:),Y(:),Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ}};
end