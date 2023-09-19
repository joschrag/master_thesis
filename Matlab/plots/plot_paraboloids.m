function [fig] = plot_paraboloids(pos_eig,a,b,V,offsets,fig,add)
%PLOT_PARABOLOIDS Summary of this function goes here
%   Detailed explanation goes here
arguments
    pos_eig (1,1) {mustBeInteger,mustBeNonnegative};
    a (1,1) {mustBeReal,mustBePositive};
    b (1,1) {mustBeReal,mustBePositive};
    V (3,3) {mustBeReal};
    offsets (3,1) {mustBeReal};
    fig matlab.ui.Figure = figure;
    add {mustBeNumericOrLogical} = true
end
if size(fig.CurrentAxes) == 0
    ax = axes(fig);
else
    ax = fig.CurrentAxes;
end

T = -5:0.1:5;
U = -5:0.1:5;
[c1,c2] = meshgrid(T,U);
if pos_eig == 1
    x =@(s,t,a,b) a.*s;
    y =@(s,t,a,b) b.*t;
    z =@(s,t,a,b) (s.^2-t.^2)./2;
else
    x =@(s,t,a,b) a.*s;
    y =@(s,t,a,b) b.*t;
    z =@(s,t,a,b)(s.^2+t.^2)./2;
end
shape = size(c1);
X = x(c1,c2,a,b);
Y = y(c1,c2,a,b);
Z = z(c1,c2,a,b);
tmp = V*[X(:),Y(:),Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
s = surf(ax,XX,YY,ZZ);
names = ["hyp Paraboloid","ell Paraboloid"];
s.DisplayName = names(pos_eig);

set(s,"EdgeColor","None")
hld = ["off","on"];
hold(ax,hld(add+1))
end

