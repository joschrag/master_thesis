function [params] = get_ellipsoid_params(a,b,c,V,offsets)
%GET_ELLIPSOID_PARAMS Summary of this function goes here
%   Detailed explanation goes here
step_arc = 0.01;
x =@(theta,phi,a,b,c) a.*sin(theta).*cos(phi);
y =@(theta,phi,a,b,c) b.*sin(theta).*sin(phi);
z =@(theta,phi,a,b,c) c.*cos(theta);
[c1,c2] = meshgrid(0:step_arc:pi+step_arc,0:step_arc:2*pi+step_arc);
shape = size(c1);
X = x(c1,c2,a,b,c);
Y = y(c1,c2,a,b,c);
Z = z(c1,c2,a,b,c);
tmp = V*[X(:),Y(:),Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
size(XX)
params = {{XX,YY,ZZ}};
end