function [params] = get_par_cylinder_params(a,b,P,offsets)
%GET_ELL_CYLINDER_PARAMS Summary of this function goes here
%   Detailed explanation goes here
U = -1:0.1:1;
V = -5:0.1:5;
[c1,c2] = meshgrid(U,V);
x =@(s,t,a) s;
y =@(s,t,a) t;
z =@(s,t,a) a.*s.^2+b.*t;
shape = size(c1);
X = x(c1,c2,a);
Y = y(c1,c2,a);
Z = z(c1,c2,a);
tmp = P*[X(:),Y(:),Z(:)]' - P*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ}};
end

