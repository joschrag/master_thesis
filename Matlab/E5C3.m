function [result] = E5C3(c)
%E3Q3_CUBIC Summary of this function goes here
%   Detailed explanation goes here
arguments       %1  2   3   4   5   6   7   8   9  10  11  12   13 14 15  16   17  18  19  20
    c (5,20) = [0, -1,  0,  0,  0,  0,  0,  0,  0,  0, -1,  6,  0, 0,  0,  0,  0, -11,  0,  7;...
                0,  0, -1,  0,  0,  0,  0,  0,  0,  0, -4,  0,  3,  0,  0,  0,  8,  0, -3, -3;...
                4,  0,  0,  0,  0,  0, 0,  -1,  0,  0,-24,  1,  0,  0,  0,  0, 44,  0,  1,-25;...
                0,  0,  0,  0,  0,  0,  0,  0, -1,  0, -2,  0,  3,  0,  0,  0, 12,  1, -4, -9;...
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  3, -1,  1, -2;...
        ];
end
% clc
c(:,10)
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
vars = [x,y,z];
var_vec = [x^3,y^3,z^3,x^2*y,x^2*z,x*y^2,x*z^2,y^2*z,y*z^2,x*y*z,x^2,y^2,z^2,x*y,x*z,y*z,x,y,z,1]';
A = -[c(:,2),c(:,3),c(:,8),c(:,9),c(:,10)*x+c(:,16)];
P = [c(:,6).*x+c(:,12),...
    c(:,7).*x+c(:,13),...
    c(:,4).*x^2+c(:,14).*x+c(:,18),...
    c(:,5).*x^2+c(:,15).*x+c(:,19),...
    c(:,1).*x^3+c(:,11).*x^2+c(:,17).*x+c(:,20)];
if rank(A) < 5
    warning("Invalid rank for matrix A!")
    return
end
P2 = inv(A)*P
lin_vars = [y^2;z^2;y;z;1]
p_var = x

[lin_vars,p_var,P2] = split_matrices(c,size(c,1),x,y,z)

cube_1 = P2(1,:)*lin_vars;   %y^3
cube_2 = P2(2,:)*lin_vars;  %z^3
quad12 = P2(3,:)*lin_vars;   %y^2*z
quad21 = P2(4,:)*lin_vars;   %y*z^2
mixed = P2(5,:)*lin_vars;    %y*z

identities = [(cube_1)*lin_vars(4) == (quad12)*lin_vars(3);...      % (y^3*)z = (y^2*z)*y
    (cube_2)*lin_vars(3) == (quad21)*lin_vars(4);...                % (z^3)*y = (z^2*y)*z
    (quad21)*lin_vars(3) == (quad12)*lin_vars(4);...                % (y^2*z)*z = (z^2*y)*y
    (quad12) == (mixed)*lin_vars(3);...                     % (y^2*z) = (y*z)*y
    (mixed)*lin_vars(4) == (quad21)];                      % (y*z)*z = (z^2*y)

sub_new = collect(identities);
old_vars = [lin_vars(3)^3,lin_vars(4)^3,lin_vars(3)^2*lin_vars(4),lin_vars(3)*lin_vars(4)^2,lin_vars(3)*lin_vars(4)];
new_vars = collect([cube_1,cube_2,quad12,quad21,mixed],[y,z]);
sub_new = subs(expand(sub_new),old_vars,new_vars);
expr = lhs(sub_new) - rhs(sub_new);
u = sym("u","real");
v = sym("v","real");

sub_vars = [y^2;z^2;y;z;1];
sub_final = subs(sub_new,sub_vars(1:2),[u;v]);
[A_,b] = equationsToMatrix(sub_final,[u;v;sub_vars(3:4)]);
A2 = [A_,-b];
vpa(simplify(A2))
coef = coeffs(det(A2),x,"All");
p_root = roots(coef);
p_root = unique(p_root(abs(imag(p_root))<10^-10));
result = [];
idx = [find(vars==p_var),find(vars==lin_vars(3)),find(vars==lin_vars(4))];
for i=1:numel(p_root)
    M = subs(A2,x,p_root(i));
    vpa(subs(det(A2),x,p_root(i)));
    if rank(M) == 5
        warning("Precision of root is too low!")
    end
    rM = rref(double(M), 1e-6);
    m = zeros(1,rank(rM));
    for j=1:rank(rM)
        m(j) = find(rM(j,:),1,"first");
    end
    r = vpa(rM(1:rank(rM),setdiff(1:5,m)));
    switch join(string(m),"")
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
            fprintf("System of equations is singular.\n");
            continue;
    end
    if ~isempty(q_root)
        for root=[q_root,r_root]'
            result=[result;[p_root(i),root']];
            result(end,:) = result(end,idx);
        end
    end
end
equations = c*var_vec;
if numel(result) > 0
    print_solutions(result,equations,x,y,z)
end
end

