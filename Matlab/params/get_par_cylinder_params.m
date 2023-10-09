function [params] = get_par_cylinder_params(a,V,offsets)
%GET_ELL_CYLINDER_PARAMS Summary of this function goes here
%   Detailed explanation goes here
T = -5:0.1:5;
U = -5:0.1:5;
[c1,c2] = meshgrid(T,U);
x =@(s,t,a) s;
y =@(s,t,a) a.*s.^2;
z =@(s,t,a) t;
shape = size(c1);
X = x(c1,c2,a);
Y = y(c1,c2,a);
Z = z(c1,c2,a);
tmp = V*[X(:),Y(:),Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ}};
end

