function [params] = get_ell_paraboloid_params(a,b,V,offsets)
%GET_ELL_PARABOLOID_PARAMS Summary of this function goes here
%   Detailed explanation goes here
T = -5:0.1:5;
U = -5:0.1:5;
[c1,c2] = meshgrid(T,U);
x =@(s,t,a,b) a.*s;
y =@(s,t,a,b) b.*t;
z =@(s,t,a,b) (s.^2-t.^2)./2;
shape = size(c1);
X = x(c1,c2,a,b);
Y = y(c1,c2,a,b);
Z = z(c1,c2,a,b);
tmp = V*[X(:),Y(:),Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ}};
end
