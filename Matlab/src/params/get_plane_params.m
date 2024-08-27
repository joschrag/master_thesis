function [params] = get_plane_params(P, offsets,opt)
    arguments
        P (3,3) {mustBeReal};
        offsets (3,1) {mustBeReal};
        opt.plotRange (1,2) {mustBeReal} = [-5,5];
    end
    T = define_plot_points(opt.plotRange);
    U = define_plot_points(opt.plotRange);
[c1,c2] = meshgrid(T,U);
y =@(s,t) s;
z =@(s,t) t;
shape = size(c1);
X = zeros(size(c1));
Y = y(c1,c2);
Z = z(c1,c2);
tmp = P*[X(:),Y(:),Z(:)]' - P*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ}};
end