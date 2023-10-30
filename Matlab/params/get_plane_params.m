function [params] = get_plane_params(V, offsets)
T = -5:0.1:5;
U = -5:0.1:5;
[c1,c2] = meshgrid(T,U);
y =@(s,t) s;
z =@(s,t) t;
shape = size(c1);
X = zeros(size(c1));
Y = y(c1,c2);
Z = z(c1,c2);
tmp = V*[X(:),Y(:),Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ}};
end