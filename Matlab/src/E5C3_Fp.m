function result = E5C3_Fp(c,prime,opt)
%E3Q3 Summary of this function goes here
%   Detailed explanation goes here
arguments                      %1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
    c (5,20) {mustBeInteger} = [1, 2, 3, 0, 0, 0, 1, 0, 0, 2, 0, 0, 1, 1, 1, 0, 1, 1, 1,-10;...
                                -1,2, 3,-1, 0, 0, 1, 0, 1, 0, 0, 3, 0, 1, 5,-1, 7, 0, 0,-10;...
                                2, 3,-1, 2, 0, 0, 1, 0, 0, 4, 1, 0, 0, 1,-3, 4,-15, 0, 0,-10;...
                                2, 2,-1, 2, 0, 3, 1, 0, 0, 0, 1, 0, 0, 1,-3, 4,-15,-3, 0,-10;...
                                2,-3,-1, 2, 0, 0, 1, 0, 4, 1, 1, 0, 0, 1,-3, 4,-15, 0, 0,-10;
        ];
    prime (1,1) {mustBeInteger} = nextprime(6);
    opt.verbose (1,1) {mustBeInteger, mustBeInRange(opt.verbose,0,2)} = 1;
    opt.error (1,1) {mustBeNumericOrLogical} = true;
    opt.log_db (1,1) {mustBeNumericOrLogical} = true;
end
t1 = tic();
x = sym("x","integer");
y = sym("y","integer");
z = sym("z","integer");
vars = [x,y,z];
all_vars = [x.^3,x.^2.*y,x.^2.*z,x.*y.^2,...
    x.*y.*z,x.*z.^2,y.^3,y.^2.*z,y.*z.^2,z.^3,x.^2,x.*y,x.*z,y.^2,y.*z,z.^2,x,y,z,1]';

[lin_vars, p_var, P2] = split_matrices_Fp(c,prime,size(c,1),x,y,z,verbose=opt.verbose,error=opt.error);

cube_1 = FF(P2.value(1,:),prime)*lin_vars;   %y^3
quad12 = FF(P2.value(2,:),prime)*lin_vars;   %y^2*z
quad21 = FF(P2.value(3,:),prime)*lin_vars;   %y*z^2
cube_2 = FF(P2.value(4,:),prime)*lin_vars;  %z^3
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
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars.value(3)),find(vars==lin_vars.value(4))]);
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
completion_time = toc(t1);
if opt.log_db
    log_to_db(c,result,completion_time,0,prime);
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

