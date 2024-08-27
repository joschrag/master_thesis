function pa_transformation(coeff_vec,opt)
%PA_TRANSFORMATION Summary of this function goes here
%   Detailed explanation goes here
%   arguments:
%       coeffs: coefficient vector of (x^2,xy,xz,y^2,yz,z^2,x,y,z,1)
arguments
    coeff_vec (1,10) {mustBeReal} = [1,1,1,1,1,1,4,-4,10,-5];
    opt.tolerance (1,1) {mustBePositive,mustBeReal} = 10^-10;
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
end
A = [coeff_vec(1),coeff_vec(2)/2,coeff_vec(3)/2;...
    coeff_vec(2)/2,coeff_vec(4),coeff_vec(5)/2;...
    coeff_vec(3)/2,coeff_vec(5)/2,coeff_vec(6)];
B = [coeff_vec(7),coeff_vec(8),coeff_vec(9)];
res = classify_wolfram(coeff_vec);
const = coeff_vec(10);

[P,lambda] = eig(A,"vector");
lambda = (abs(lambda) > opt.tolerance).*lambda;
P = [P(:,1)./norm(P(:,1)),P(:,2)./norm(P(:,2)),P(:,3)./norm(P(:,3))];
d_vec = zeros(1,2);
d_ind = 1;
offsets = zeros(3,1);
if any(B~=0)
    lin_coeffs = B*P;
    for i=1:3
        if lambda(i) ~= 0
            [offsets(i),const] = square_complete(lambda(i),lin_coeffs(i),const);
        else
            d_vec(d_ind) = lin_coeffs(i);
            d_ind = d_ind + 1;
        end
    end
end
d0 = d_vec(1);
d1 = d_vec(2);
switch res
    case 3
        params = ellipsoid_nf(lambda,P,offsets,const);
        dname = "ellipsoid";
    case 4
        params = hyperboloid1_nf(lambda,P,offsets,const);
        dname = "hyperboloid1";
    case 5
        params = hyperboloid2_nf(lambda,P,offsets,const);
        dname = "hyperboloid2";
    case 6
        params = ell_paraboloid_nf(lambda,P,offsets,d0,const);
        dname = "ell paraboloid";
    case 7
        params = hyp_paraboloid_nf(lambda,P,offsets,d0,const);
        dname = "hyp paraboloid";
    case 2
        params = ell_cone_nf(lambda,P,offsets);
        dname = "ell cone";
    case 8
        params = ell_cylinder_nf(lambda,P,offsets,const);
        dname = "ell cylinder";
    case 10
        params = hyp_cylinder_nf(lambda,P,offsets,const);
        dname = "hyp cylinder";
    case 11
        params = par_cylinder_nf(lambda,P,offsets,d0,d1,const);
        dname = "par cylinder";
    case 12
        params = par_planes_nf(lambda,P,offsets,const);
        dname = "par planes";
    case 13
        params = plane_nf(lambda,P,offsets);
        dname = "plane";
    case 14
        params = cross_planes_nf(lambda,P,offsets);
        dname = "cross planes";
    case 9
        params = line_nf(lambda,P,offsets);
        dname = "line";
    case 1
        params = {{0,0,0}};
        dname = "point";
end

ax = gca();
for i = 1:numel(params)
    coords = params{i};
    if all(size(coords{3})>1)
        go = surf(ax,coords{1},coords{2},coords{3});
        set(go,"FaceColor","b")
        set(go,"FaceAlpha",0.4)
        set(go,"EdgeColor","None")
    elseif any(size(coords{3})>1)
        go = plot3(ax,coords{1},coords{2},coords{3});
    else
        go = scatter3(ax,coords{1},coords{2},coords{3});
    end
    go.DisplayName = dname;
    hold(ax,"on")
end
xlabel("X")
ylabel("Y")
zlabel("Z")
end