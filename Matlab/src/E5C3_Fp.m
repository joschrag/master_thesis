function result = E5C3_Fp(c,prime)
%E3Q3 Summary of this function goes here
%   Detailed explanation goes here
arguments                      %1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
    c (5,20) {mustBeInteger} = [1, 2, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1,-10;...
                                -1,2, 3,-1, 0, 0, 1, 0, 1, 0, 0, 3, 0, 1, 5,-1, 7, 0, 0,-10;...
                                2, 3,-1, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1,-3, 4,-15, 0, 0,-10;...
                                2, 2,-1, 2, 2, 3, 1, 0, 0, 0, 1, 0, 0, 1,-3, 4,-15,-3, 0,-10;...
                                2,-3,-1, 2, 0, 0, 1, 0, 4, 0, 1, 0, 0, 1,-3, 4,-15, 0, 0,-10;
        ];
    prime (1,1) {mustBeInteger} = nextprime(6);
end
x = sym("x","integer");
y = sym("y","integer");
z = sym("z","integer");
vars = [x,y,z];
all_vars = [x^3,y^3,z^3,x^2*y,x^2*z,x*y^2,y^2*z,x*z^2,y*z^2,x*y*z,x^2,y^2,z^2,x*y,x*z,y*z,x,y,z,1]';
A = FF(-[c(:,2),c(:,3),c(:,7),c(:,9),c(:,10)*x+c(:,16)],prime);
[s1,~] = size(c);
r_A = rank(A);

if r_A == 5
    fprintf("Using P(x)\n")
    Q = A;
    P = FF([c(:,6).*x+c(:,12),...
        c(:,8).*x+c(:,13),...
        c(:,4).*x^2+c(:,14).*x+c(:,18),...
        c(:,5).*x^2+c(:,15).*x+c(:,19),...
        c(:,1).*x^3+c(:,11).*x^2+c(:,17).*x+c(:,20)],prime);
    lin_vars = FF([y^2;z^2;y;z;1],prime);
    p_var = x;
else
    B = -FF([c(:,1),c(:,3),c(:,5),c(:,8),c(:,10)*y+c(:,15)],prime);
    r_B = rank(B);
    if r_B == 5
        fprintf("Using P(y)\n")
        Q = B;
        P = FF([c(:,4).*y+c(:,14),...
            c(:,9).*y+c(:,16),...
            c(:,6).*y^2+c(:,14).*y+c(:,17),...
            c(:,7).*y^2+c(:,16).*y+c(:,19),...
            c(:,2).*y^3+c(:,12).*y^2+c(:,18).*y+c(:,20)],prime);
        lin_vars = FF([x^2;z^2;x;z;1],prime);
        p_var = y;
    else
        C = -FF([c(:,1),c(:,2),c(:,4),c(:,6),c(:,10)*z+c(:,14)],prime);
        r_C = rank(C);
        if r_C == 5
            fprintf("Using P(z)\n")
            Q = C;
            P = FF([c(:,5).*z+c(:,11),...
                c(:,7).*z+c(:,12),...
                c(:,8).*z^2+c(:,15).*z+c(:,17),...
                c(:,9).*z^2+c(:,16).*z+c(:,18),...
                c(:,3).*z^3+c(:,13).*z^2+c(:,19).*z+c(:,20)],prime);
            lin_vars = FF([x^2;y^2;x;y;1],prime);
            p_var = z;
        else
            %warning("No valid matrix found!")
            result = [];
            return;
        end
    end
end

P2 = inv(Q)*P;
% fprintf("P2:\n%s %s %s %s %s\n%s %s %s %s %s\n%s %s %s %s %s\n%s %s %s %s %s\n%s %s %s %s %s\n",string(P2.value))
cube_1 = FF(P2.value(1,:),prime)*lin_vars;   %y^3
cube_2 = FF(P2.value(2,:),prime)*lin_vars;  %z^3
quad12 = FF(P2.value(3,:),prime)*lin_vars;   %y^2*z
quad21 = FF(P2.value(4,:),prime)*lin_vars;   %y*z^2
mixed = FF(P2.value(5,:),prime)*lin_vars;    %y*z

identities = [(cube_1)*FF(lin_vars.value(4),prime) - (quad12)*FF(lin_vars.value(3),prime);...      % (y^3*)z = (y^2*z)*y
    (cube_2)*FF(lin_vars.value(3),prime) - (quad21)*FF(lin_vars.value(4),prime);...                % (z^3)*y = (z^2*y)*z
    (quad21)*FF(lin_vars.value(3),prime) - (quad12)*FF(lin_vars.value(4),prime);...                % (y^2*z)*z = (z^2*y)*y
    (quad12) - (mixed)*FF(lin_vars.value(3),prime);...                     % (y^2*z) = (y*z)*y
    (mixed)*FF(lin_vars.value(4),prime) - (quad21)];                      % (y*z)*z = (z^2*y)



sub1 = FF.empty(5,0);
sub2 = FF.empty(5,0);
sub3 = FF.empty(5,0);
for i = 1:5
    sub1(i) = collect(identities(i));
end
old_vars = [lin_vars.value(3)^3,...
    lin_vars.value(4)^3,...
    lin_vars.value(3)^2*lin_vars.value(4),...
    lin_vars.value(3)*lin_vars.value(4)^2,...
    lin_vars.value(3)*lin_vars.value(4)];
new_vars = [cube_1.value,cube_2.value,quad12.value,quad21.value,mixed.value];
for i = 1:5
    sub2(i) = subs(sub1(i),old_vars,new_vars);
end
for i = 1:5
    ex = expand(sub2(i));
    sub3(i) = subs(ex,old_vars,new_vars);
end
eq = [sub3(1).value,sub3(2).value,sub3(3).value,sub3(4).value,sub3(5).value];
eq = simplify(eq);
u = sym("u","real");
v = sym("v","real");

sub_final = subs(eq,lin_vars.value(1:2),[u;v]);
[A_,b] = equationsToMatrix(sub_final,[u;v;lin_vars.value(3:4)]);
A2 = FF([A_,-b],prime);
% fprintf("A2:\n%s %s %s %s %s\n%s %s %s %s %s\n%s %s %s %s %s\n%s %s %s %s %s\n%s %s %s %s %s\n",string(A2.value))
pol = coeffs(det(A2),p_var,"All");
p_root = get_gf_root(pol,prime);
if numel(p_root) == 0
    fprintf("No solutions to the equations.\n")
    result = [];
    return;
end
idx = [find(vars==p_var),find(vars==lin_vars.value(3)),find(vars==lin_vars.value(4))];
result = [];
for i = 1:numel(p_root)
    p_0 = p_root(i);
    M = FF(vpa(subs(A2.value,p_var,p_0)),prime);
    % fprintf("p_0: %i\n",p_0)
    % fprintf("M:\n%i %i %i %i %i\n%i %i %i %i %i\n%i %i %i %i %i\n%i %i %i %i %i\n%i %i %i %i %i\n\n",M.value)
    if M.value == zeros(5)
        fprintf("Matrix is singular.\n");
        continue;
    end
    rM = rref(M);
    m = zeros(1,rank(rM));
    for j=1:rank(rM)
        m(j) = find(rM(j,:),1,"first");
    end
    r = vpa(rM(1:rank(rM),setdiff(1:5,m)));
    switch join(string(m),"")
        case "1"
            [q_root,r_root] = rank1_1_fp(r,prime);
        case "2"
            [q_root,r_root] = rank1_2_fp(r,prime);
        case "3"
            [q_root,r_root] = rank1_3_fp(r,prime);
        case "4"
            [q_root,r_root] = rank1_4_fp(r,prime);
        case "12"
            [q_root,r_root] = rank2_12_fp(r,prime);
        case "13"
            [q_root,r_root] = rank2_13_fp(r,prime);
        case "14"
            [q_root,r_root] = rank2_14_fp(r,prime);
        case "23"
            [q_root,r_root] = rank2_23_fp(r,prime);
        case "24"
            [q_root,r_root] = rank2_24_fp(r,prime);
        case "34"
            [q_root,r_root] = rank2_34_fp(r,prime);
        case "123"
            [q_root,r_root] = rank3_123_fp(r,prime);
        case "124"
            [q_root,r_root] = rank3_124_fp(r,prime);
        case "134"
            [q_root,r_root] = rank3_134_fp(r,prime);
        case "234"
            [q_root,r_root] = rank3_234_fp(r,prime);
        case "1234"
            [q_root,r_root] = rank4_1234_fp(r,prime);
        otherwise
            fprintf("System of equations is singular.\n");
            continue;
    end
    if ~isempty(q_root)
        for root=[q_root,r_root]'
            result=[result;[p_0,root']];
            result(end,:) = result(end,idx);
        end
    end
end
if isempty(result)
    result = [];
    return
end
equations = c*all_vars;
if numel(result) > 0
    print_solutions(result,equations,x,y,z,prime)
end
end

