function plot_and_color_implicit(result, m, c, opt)
arguments
    result (:,3) {mustBeReal};
    m (1,1) {mustBeInteger,mustBePositive}
    c (:,20) {mustBeReal};
    opt.intervals = [-5,5];
    opt.density (1,5) = ones(1,5).*100;
end
color_list = ["#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999"];
if ~isempty(result)
    min_x = min(result(:,1));
    max_x = max(result(:,1));
    min_y = min(result(:,2));
    max_y = max(result(:,2));
    min_z = min(result(:,3));
    max_z = max(result(:,3));
    min_dist = 10;
    extra = abs(max(max_x-min_x,min_dist))*0.5;
    ex_x_roots = [min_x-extra,max_x+extra];
    extra = abs(max(max_y-min_y,min_dist))*0.5;
    ex_y_roots = [min_y-extra,max_y+extra];
    extra = abs(max(max_z-min_z,min_dist))*0.5;
    ex_z_roots = [min_z-extra,max_z+extra];
    interval = [ex_x_roots,ex_y_roots,ex_z_roots];
else
    interval = opt.intervals;
end
figure
for i=1:m
    f = @(x,y,z) sum(repmat(double(c(i,:))',size(x)).*[x.^3; x.^2.*y; x.^2.*z; ...
                                         x.*y.^2; x.*y.*z; x.*z.^2; ...
                                         y.^3; y.^2.*z; y.*z.^2; ...
                                         z.^3; x.^2; x.*y; x.*z; ...
                                         y.^2; y.*z; z.^2; x; y; z; ones(size(x))], 1);
    s = fimplicit3(f,double(interval),"MeshDensity",opt.density(i));
    set(s,"FaceAlpha",0.4,"FaceColor",color_list(i),"EdgeColor","none");
    hold on
end
for i=1:size(result,1)
    scatter3(result(i,1),result(i,2),result(i,3),50,"black","filled")
end
xlabel("x")
ylabel("y")
zlabel("z")
end