function params = get_line_params(P,offsets,opt)
    arguments
        P (3,3) {mustBeReal};
        offsets (3,1) {mustBeReal};
        opt.plotRange (1,2) {mustBeReal} = [-5,5];
    end
    U = define_plot_points(opt.plotRange);
shape = size(U);
X = zeros(shape);
Y = zeros(shape);
Z = U;
tmp = P*[X(:),Y(:),Z(:)]' - P*offsets;
XX = reshape(tmp(1,:),shape);
YY = reshape(tmp(2,:),shape);
ZZ = reshape(tmp(3,:),shape);
params = {{XX,YY,ZZ}};
end