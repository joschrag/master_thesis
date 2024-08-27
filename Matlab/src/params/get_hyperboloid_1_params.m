function [params] = get_hyperboloid_1_params(a,b,c,P,offsets,opt)
%GET_HYPERBOLOID_1_PARAMS Summary of this function goes here
%   Detailed explanation goes here
arguments
    a (1,1) {mustBeReal};
    b (1,1) {mustBeReal};
    c (1,1) {mustBeReal};
    P (3,3) {mustBeReal};
    offsets (3,1) {mustBeReal};
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
end
step_arc = 0.01;
x =@(s,t,a,b,c)  a.*cosh(s).*cos(t);
y =@(s,t,a,b,c) b.*cosh(s).*sin(t);
z =@(s,t,a,b,c) c.*sinh(s);
U = define_plot_points(opt.plotRange);
[c1,c2] = meshgrid(U,0:step_arc:2*pi+step_arc);
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