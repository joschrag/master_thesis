function [params] = get_ell_paraboloid_params(a,b,V,offsets,dir)
%GET_ELL_PARABOLOID_PARAMS Summary of this function goes here
%   Detailed explanation goes here
T = [0:0.01:0.4,0.5:0.1:5];
step_arc = 0.01;
U = -pi:step_arc:pi+step_arc;
[c1,c2] = meshgrid(T,U);
x =@(s,t,a,b) a.*s.*cos(t);
y =@(s,t,a,b) b.*s.*sin(t);
z =@(s,t,a,b) s.^2;
shape = size(c1);
X = x(c1,c2,a,b);
Y = y(c1,c2,a,b);
Z = dir.*z(c1,c2,a,b);
tmp = V*[X(:),Y(:),Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ}};
end

