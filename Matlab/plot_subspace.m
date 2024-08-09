function [] = plot_subspace(u_root,v_root,r,base)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
    u_root (:,:) {mustBeReal};
    v_root (:,:) {mustBeReal};
    r (1,1) {mustBeInteger,mustBePositive};
    base (:,:) {mustBeReal};
end
figure
for i=1:size(u_root,2)
    % cv = mean(v_root(:,i));
    % cu = mean(u_root(:,i));
    % [~, order] = sort(atan2(v_root(:,i)-cv, u_root(:,i)-cu));
    % u_root(:,i) = u_root(order,i);
    % v_root(:,i) = v_root(order,i);
    scatter3(u_root(:,i),v_root(:,i),u_root(:,i).^2,25,'MarkerEdgeColor',"#D73027",...
        'MarkerFaceColor',"#D73027")
    hold on
    scatter3(u_root(:,i),v_root(:,i),v_root(:,i).^2,25,'MarkerEdgeColor',"#4575B4",...
        'MarkerFaceColor',"#4575B4")
end
if ~isempty(u_root) && ~isempty(v_root)
    min_u = min(u_root,[],"all");
    max_u = max(u_root,[],"all");
    min_v = min(v_root,[],"all");
    max_v = max(v_root,[],"all");
    min_dist = 5;
    extra = abs(max(max_u-min_u,min_dist))*0.35;
    ex_u_roots = [min_u-extra,max_u+extra];
    extra = abs(max(max_v-min_v,min_dist))*0.35;
    ex_v_roots = [min_v-extra,max_v+extra];
else
    ex_u_roots = [-5,5];
    ex_v_roots = [-5,5];
end
[U,V] = meshgrid(ex_u_roots(1):0.1:ex_u_roots(2),ex_v_roots(1):0.1:ex_v_roots(2));
U = double(U);
V = double(V);
if r == 1
    Z =@(i) base(i,1).*V.^2+base(i,2).*U+base(i,3).*V+base(i,4);
elseif r == 2
    Z =@(i) base(i,1).*U+base(i,2).*V+base(i,3);
elseif r == 3
    Z =@(i) base(i,1).*V+base(i,2);
else
    Z =@(i) 0.*V+base(i,1);
end
s = surf(U,V,double(Z(1)),"EdgeColor","none");
set(s,'FaceAlpha',0.8,'FaceColor',"#FEE090",'EdgeColor',"none");
hold on
s = surf(U,V,U.^2,"EdgeColor","none");
set(s,'FaceAlpha',0.8,'FaceColor',"#FC8D59",'EdgeColor',"none");
s = surf(U,V,double(Z(2)),"EdgeColor","none");
set(s,'FaceAlpha',0.8,'FaceColor',"#91BFDB",'EdgeColor',"none");
hold on
s = surf(U,V,V.^2,"EdgeColor","none");
set(s,'FaceAlpha',0.8,'FaceColor',"#E0F3F8",'EdgeColor',"none");
axis(double([ex_u_roots,ex_v_roots]))
if r~= 1
    for i=1:numel(u_root)
        plot3([u_root(i),u_root(i),],[v_root(i),v_root(i)],zlim,"--k");
    end
end
xlabel("u")
ylabel("v")
end