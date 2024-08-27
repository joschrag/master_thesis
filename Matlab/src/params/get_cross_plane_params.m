function [params] = get_cross_plane_params(a,b, P, offsets, opt)
%GET_PARALLEL_PLANE_PARAMS Summary of this function goes here
%   Detailed explanation goes here
arguments
    a (1,1) {mustBeReal};
    b (1,1) {mustBeReal};
    P (3,3) {mustBeReal};
    offsets (3,1) {mustBeReal};
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
end
T = define_plot_points(opt.plotRange);
U = define_plot_points(opt.plotRange);
[c1,c2] = meshgrid(T,U);
shape = size(c1);
x =@(s,t,a) s;
y =@(s,t,a) b./a.*s;
z =@(s,t,a) t;
X = x(c1,c2,a);
Y = y(c1,c2,a);
Z = z(c1,c2,a);
tmp = P*[X(:),Y(:),Z(:)]' - P*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
tmp = P*[-X(:),Y(:),Z(:)]' - P*offsets;
XX2 = reshape(tmp(1,:),shape);
YY2 = reshape(tmp(2,:),shape);
ZZ2 = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ},{XX2,YY2,ZZ2}};
end


