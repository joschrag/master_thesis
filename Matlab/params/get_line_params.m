function params = get_line_params(V,offsets)
U = -5:0.1:5;
shape = size(U);
X = zeros(shape);
Y = zeros(shape);
Z = U;
tmp = V*[X(:),Y(:),Z(:)]' - V*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ}};
end