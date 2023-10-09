function fig = plot_ell_cone(a,b,c,V,offsets,fig,add)
arguments
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
s = 0:0.1:2*pi+0.1;
t = 0:10;
x =@(s,t,a,b,c) a.*t.*cos(s);
y =@(s,t,a,b,c) b.*t.*sin(s);
z =@(s,t,a,b,c) c.*t;
[p1,p2] = meshgrid(s,t);
X = x(p1,p2,a,b,c);
Y = y(p1,p2,a,b,c);
Z = z(p1,p2,a,b,c);
shape = size(p1);
tmp = V*[X(:),Y(:),Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
s = surf(ax,XX,YY,ZZ);
s.EdgeColor = 'none';
s.DisplayName = "Ell Cone";
hold(ax,"on")
tmp = V*[X(:),Y(:),-Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
s = surf(ax,XX,YY,ZZ);
s.EdgeColor = 'none';
s.DisplayName = "Ell Cone";
hld = ["off","on"];
hold(ax,hld(add+1))
end