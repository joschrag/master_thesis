function [params] = get_ell_cone_params(a,b,c,P,offsets,opt)
%GET_ELL_CONE_PARAMS Summary of this function goes here
%   Detailed explanation goes here
arguments
    a (1,1) {mustBeReal};
    b (1,1) {mustBeReal};
    c (1,1) {mustBeReal};
    P (3,3) {mustBeReal};
    offsets (3,1) {mustBeReal};
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
end
t = define_plot_points([0,opt.plotRange(2)]);
s = 0:0.1:2*pi+0.1;
x =@(s,t,a,b,c) a.*t.*cos(s);
y =@(s,t,a,b,c) b.*t.*sin(s);
z =@(s,t,a,b,c) c.*t;
[p1,p2] = meshgrid(s,t);
X = x(p1,p2,a,b,c);
Y = y(p1,p2,a,b,c);
Z = z(p1,p2,a,b,c);
shape = size(p1);
tmp = P*[X(:),Y(:),Z(:)]' - P*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
tmp = P*[X(:),Y(:),-Z(:)]' - P*offsets;
XX2 = reshape(tmp(1,:),shape);
YY2 = reshape(tmp(2,:),shape);
ZZ2 = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ},{XX2,YY2,ZZ2}};
end

