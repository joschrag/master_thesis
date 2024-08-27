function [params] = get_ell_paraboloid_params(a,b,P,offsets,dir,opt)
%GET_ELL_PARABOLOID_PARAMS Summary of this function goes here
%   Detailed explanation goes here
arguments
    a (1,1) {mustBeReal};
    b (1,1) {mustBeReal};
    P (3,3) {mustBeReal};
    offsets (3,1) {mustBeReal};
    dir (1,1) {mustBeMember(dir,[-1,1])};
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
end
T = define_plot_points(opt.plotRange);
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
tmp = P*[X(:),Y(:),Z(:)]' - P*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ}};
end

