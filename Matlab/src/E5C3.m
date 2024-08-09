function [result] = E5C3(c)
%E3Q3_CUBIC Summary of this function goes here
%   Detailed explanation goes here
arguments       %1  2   3   4   5   6   7   8   9  10  11  12   13 14 15  16   17  18  19  20
    c (5,20) = [1,-3,-2,0,0,0,-1,4,-1,-2,2,6,2,-1,-4,5,-1,-5,2,-2;...
        1,2,-1,-3,0,-1,0,1,-2,2,-4,2,4,3,4,-5,2,-8,2,1;...
        -2,2,1,-4,0,2,-5,0,2,2,1,2,-8,22,-6,-8,9,-20,13,-3;...
        -2,2,-2,-1,0,2,2,-2,-2,2,6,-2,1,-4,9,-5,-4,-4,-5,9;...
        0,-2,0,-1,0,0,-3,0,4,0,-2,6,1,6,-8,-2,3,2,3,-7;...
        ];
end
% close("all")
% clearvars -except c
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
vars = [x,y,z];
var_vec = [x.^3,x.^2.*y,x.^2.*z,x.*y.^2,...
    x.*y.*z,x.*z.^2,y.^3,y.^2.*z,y.*z.^2,z.^3,x.^2,x.*y,x.*z,y.^2,y.*z,z.^2,x,y,z,1]';
m = size(c,1);
[lin_vars,p_var,P2] = split_matrices(c,m,x,y,z);

cube_1 = P2(1,:)*lin_vars;   %y^3
quad12 = P2(2,:)*lin_vars;   %y^2*z
quad21 = P2(3,:)*lin_vars;   %y*z^2
cube_2 = P2(4,:)*lin_vars;  %z^3
mixed = P2(5,:)*lin_vars;    %y*z

identities = [(cube_1)*lin_vars(4) - (quad12)*lin_vars(3);...      % (y^3*)z = (y^2*z)*y
    (cube_2)*lin_vars(3) - (quad21)*lin_vars(4);...                % (z^3)*y = (z^2*y)*z
    (quad21)*lin_vars(3) - (quad12)*lin_vars(4);...                % (y^2*z)*z = (z^2*y)*y
    (quad12) - (mixed)*lin_vars(3);...                     % (y^2*z) = (y*z)*y
    (mixed)*lin_vars(4) - (quad21)];                      % (y*z)*z = (z^2*y)

sub_new = collect(identities);
old_vars = [lin_vars(3)^3,lin_vars(4)^3,lin_vars(3)^2*lin_vars(4),lin_vars(3)*lin_vars(4)^2,lin_vars(3)*lin_vars(4)];
new_vars = collect([cube_1,cube_2,quad12,quad21,mixed],lin_vars(3:4));
sub_new = expand(subs(expand(sub_new),old_vars,new_vars));
u = sym("u","real");
v = sym("v","real");
sub_final = subs(sub_new,lin_vars(1:2),[u;v]);
[A_,b] = equationsToMatrix(sub_final,[u;v;lin_vars(3:4)]);
A2 = [A_,-b];
coef = coeffs(det(A2),p_var,"All");
p_root = roots(coef);
p_root = unique(p_root(abs(imag(p_root))<10^-10));
result = [];
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars(3)),find(vars==lin_vars(4))]);
for i=1:numel(p_root)
    M = subs(A2,p_var,p_root(i));
    % vpa(subs(det(A2),p_var,p_root(i)))
    if rank(M) == 5
        warning("Precision of root is too low!")
    end
    rM = rref(double(M), 1e-6);
    col = zeros(1,rank(rM));
    for j=1:rank(rM)
        col(j) = find(rM(j,:),1,"first");
    end
    r = vpa(rM(1:rank(rM),setdiff(1:5,col)));
    p_root(i)
    col

    switch join(string(col),"")
        case "1"
            [q_root,r_root] = rank1_1(r);
        case "2"
            [q_root,r_root] = rank1_2(r);
        case "3"
            [q_root,r_root] = rank1_3(r);
        case "4"
            [q_root,r_root] = rank1_4(r);
        case "12"
            [q_root,r_root] = rank2_12(r);
        case "13"
            [q_root,r_root] = rank2_13(r);
        case "14"
            [q_root,r_root] = rank2_14(r);
        case "23"
            [q_root,r_root] = rank2_23(r);
        case "24"
            [q_root,r_root] = rank2_24(r);
        case "34"
            [q_root,r_root] = rank2_34(r);
        case "123"
            [q_root,r_root] = rank3_123(r);
        case "124"
            [q_root,r_root] = rank3_124(r);
        case "134"
            [q_root,r_root] = rank3_134(r);
        case "234"
            [q_root,r_root] = rank3_234(r);
        case "1234"
            [q_root,r_root] = rank4_1234(r);
        otherwise
            % fprintf("System of equations is singular.\n");
            continue;
    end
    if ~isempty(q_root)
        for root=[q_root,r_root]'
            result=[result;[p_root(i),root']];
            result(end,:) = result(end,idx);
        end
        % plot_subspace(M)
    end
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
    interval = [-5,5];
end
dens = [25,25,25,25,25];
figure
for i=1:m
    f =@(u,v,w) double(subs(c(i,:)*var_vec,{x,y,z},{u,v,w}));
    s = fimplicit3(f,double(interval),"MeshDensity",dens(i));
    set(s,"FaceAlpha",0.4,"FaceColor",color_list(i),"EdgeColor","none");
    hold on
end
for i=1:size(result,1)
    scatter3(result(i,1),result(i,2),result(i,3),50,"black","filled")
end
xlabel("x")
ylabel("y")
zlabel("z")
equations = c*var_vec;
if numel(result) > 0
    print_solutions(result,equations,x,y,z)
end
end

