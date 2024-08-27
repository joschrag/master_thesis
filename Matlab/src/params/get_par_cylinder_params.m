function [params] = get_par_cylinder_params(a,b,P,offsets,opt)
%GET_ELL_CYLINDER_PARAMS Summary of this function goes here
%   Detailed explanation goes here
arguments
    a (1,1) {mustBeReal};
    b (1,1) {mustBeReal};
    P (3,3) {mustBeReal};
    offsets (3,1) {mustBeReal};
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
end
T = -1:0.1:1;
U = define_plot_points(opt.plotRange);
[c1,c2] = meshgrid(T,U);
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

