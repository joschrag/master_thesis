function fig = plot_quadric_3d(pos_eig,neg_eig,a,b,c,V,offsets,fig,add)
arguments
    pos_eig (1,1) {mustBeInteger,mustBeNonnegative};
    neg_eig (1,1) {mustBeInteger,mustBeNonnegative};
    a (1,1) {mustBeReal,mustBePositive};
    b (1,1) {mustBeReal,mustBePositive};
    c (1,1) {mustBeReal,mustBePositive};
    V (3,3) {mustBeReal};
    offsets (3,1) {mustBeReal};
    fig matlab.ui.Figure = gcf;
    add {mustBeNumericOrLogical} = true
end
if size(fig.CurrentAxes) == 0
    ax = axes(fig);
else
    ax = fig.CurrentAxes;
end
step_arc = 0.01;
if pos_eig == 3 % Ellipsoid
    x =@(theta,phi,a,b,c) a.*sin(theta).*cos(phi);
    y =@(theta,phi,a,b,c) b.*sin(theta).*sin(phi);
    z =@(theta,phi,a,b,c) c.*cos(theta);
    [c1,c2] = meshgrid(0:step_arc:pi+step_arc,0:step_arc:2*pi+step_arc);
elseif pos_eig == 2 && neg_eig == 1 % 1 schalig Hyperboloid
    x =@(s,t,a,b,c)  a.*cosh(s).*cos(t);
    y =@(s,t,a,b,c) b.*cosh(s).*sin(t);
    z =@(s,t,a,b,c) c.*sinh(s);
    [c1,c2] = meshgrid(0:step_arc:2*pi+step_arc,-10:0.1:10);
elseif pos_eig == 1 && neg_eig == 2 % 2 schalig Hyperboloid
    x =@(s,t,a,b,c) a.*sinh(s).*cos(t);
    y =@(s,t,a,b,c) b.*sinh(s).*sin(t);
    z =@(s,t,a,b,c) c.*cosh(s);
    [c1,c2] = meshgrid(0:step_arc:2*pi+step_arc,-10:0.1:10);
end

shape = size(c1);
X = x(c1,c2,a,b,c);
Y = y(c1,c2,a,b,c);
Z = z(c1,c2,a,b,c);
tmp = V*[X(:),Y(:),Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
s = surf(ax,XX,YY,ZZ);
set(s,"EdgeColor","None")
if pos_eig == 3 % Ellipsoid
    s.DisplayName = "Ellipsoid";
elseif pos_eig == 2 && neg_eig == 1 % 1 schalig Hyperboloid
    s.DisplayName = "1 Hyperboloid";
elseif pos_eig == 1 && neg_eig == 2 % 2 schalig Hyperboloid
    s.DisplayName = "2 Hyperboloid";
end


if pos_eig == 1 && neg_eig == 2
    hold(ax,"on")
    tmp = V*[X(:),Y(:),-Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
    s = surf(ax,XX,YY,ZZ);
    set(s,"EdgeColor","None")
    s.DisplayName = "2 Hyperboloid";
end
hld = ["off","on"];
hold(ax,hld(add+1))
end