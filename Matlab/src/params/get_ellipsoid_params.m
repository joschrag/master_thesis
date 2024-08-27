function [params] = get_ellipsoid_params(a,b,c,P,offsets)
%GET_ELLIPSOID_PARAMS Summary of this function goes here
%   Detailed explanation goes here
arguments
    a (1,1) {mustBeReal};
    b (1,1) {mustBeReal};
    c (1,1) {mustBeReal};
    P (3,3) {mustBeReal};
    offsets (3,1) {mustBeReal};
end
step_arc = 0.01;
x =@(theta,phi,a,b,c) a.*sin(theta).*cos(phi);
y =@(theta,phi,a,b,c) b.*sin(theta).*sin(phi);
z =@(theta,phi,a,b,c) c.*cos(theta);
[c1,c2] = meshgrid(0:step_arc:pi+step_arc,0:step_arc:2*pi+step_arc);
shape = size(c1);
X = x(c1,c2,a,b,c);
Y = y(c1,c2,a,b,c);
Z = z(c1,c2,a,b,c);
tmp = P*[X(:),Y(:),Z(:)]' - P*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ}};
end