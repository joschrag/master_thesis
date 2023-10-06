function fig = plot_parabola(pos_eig,a,V,offsets,fig,add)
%PLOT_PARABOLA Summary of this function goes here
%   Detailed explanation goes here
arguments
    pos_eig (1,1) {mustBeInteger,mustBeNonnegative};
    a (1,1) {mustBeReal,mustBePositive};
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
x =@(s,t,a) s;
y =@(s,t,a) a.*s.^2;
z =@(s,t,a) t;
shape = size(c1);
X = x(c1,c2,a);
Y = y(c1,c2,a);
Z = z(c1,c2,a);
tmp = V*[X(:),Y(:),Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
s = surf(ax,XX,YY,ZZ);
s.DisplayName = "Parabola";
set(s,"EdgeColor","None")
hld = ["off","on"];
hold(ax,hld(add+1))
end

