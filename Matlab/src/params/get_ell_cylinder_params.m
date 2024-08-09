function [params] = get_ell_cylinder_params(a,b,P,offsets)
%GET_ELL_CYLINDER_PARAMS Summary of this function goes here
%   Detailed explanation goes here
U = -5:0.1:5;
step_arc = 0.01;
V = 0:step_arc:2*pi+step_arc;
[c1,c2] = meshgrid(U,V);
x =@(s,t,a,b) a.*cos(t);
y =@(s,t,a,b) b.*sin(t);
z =@(s,t,a,b) s;
shape = size(c1);
X = x(c1,c2,a,b);
Y = y(c1,c2,a,b);
Z = z(c1,c2,a,b);
tmp = P*[X(:),Y(:),Z(:)]' - P*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ}};
end