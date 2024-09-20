function [result] = E4C3_Fp(c,prime,opt)
arguments       %1  2   3   4   5   6   7   8   9  10  11  12   13 14 15  16   17  18  19  20
    c (4,20) = [0,-3,-2,0,0,0,-1,0,0,0,0,6,2,0,0,0,-15,-5,2,-5;...
        0,2,-1,0,0,0,0,-1,0,0,0,2,4,0,0,0,10,-8,2,-2;...
        -1,2,1,0,0,0,0,0,-1,0,0,2,-8,0,0,0,3,7,13,-4;...
        1,2,-2,0,0,0,0,0,0,-1,0,-2,1,0,0,0,-7,-4,-5,-5;...
        ];
    prime (1,1) {mustBeInteger} = nextprime(6);
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 1;
    opt.error (1,1) {mustBeNumericOrLogical} = true;
    opt.log_db (1,1) {mustBeNumericOrLogical} = true;
end
t1 = tic();
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
vars = [x,y,z];
var_vec = [x.^3,x.^2.*y,x.^2.*z,x.*y.^2,...
    x.*y.*z,x.*z.^2,y.^3,y.^2.*z,y.*z.^2,z.^3,x.^2,x.*y,x.*z,y.^2,y.*z,z.^2,x,y,z,1]';
m = size(c,1);
complete_idx = {[4,14],[5,15],[6,16],[2,12,18],[3,13,19],[1,11,17,20];...
    [2,11],[5,13],[9,16],[4,12,17],[8,15,19],[7,14,18,20];...
    [3,11],[5,12],[8,14],[6,13,17],[9,15,18],[10,16,19,20]};
q_idx = [7,8,9,10;1,3,6,10;1,2,4,7];    %indices of cubic monomials of lin_vars
Q = cell(1,3);
ranks = zeros(1,3);
for i=1:3
    Q{i} = FF(-[c(:,q_idx(i,1)),c(:,q_idx(i,2)),c(:,q_idx(i,3)),c(:,q_idx(i,4))],prime);
    ranks(i) = rank(Q{i});
end
result = [];

if any(ranks==4)
    I = find(ranks==4);
else
    return
end

valid_P2 = false;
if ~isempty(I)
    for idx = I
        var_idx = setdiff(1:3,idx);
        all_lin_vars = [vars(var_idx(1)).^2,vars(var_idx(1))*vars(var_idx(2)),vars(var_idx(2)).^2,vars(setdiff(1:3,I)),1]';
        if opt.verbose > 0
            fprintf("Using P(%s)\n",string(vars(idx)))
        end
        p_var = vars(idx);
        p_var_pow = [1;p_var;p_var^2;p_var^3];
        P = sym.zeros(4,6);
        for j=1:6
            P(:,j) = c(:,complete_idx{idx,j})*p_var_pow(numel(complete_idx{idx,j}):-1:1);
        end
        lin_vars = all_lin_vars;
        P2 = Q{idx}^(-1)*FF(P,prime);
        if all(all(P2.value(:,1:1:3)==sym(0)))
            valid_P2 = true;
            break
        end
    end
else
    return
end

if ~valid_P2
    return
end


cube_1 = FF(P2.value(1,:),prime)*FF(lin_vars,prime);   %y^3
quad12 = FF(P2.value(2,:),prime)*FF(lin_vars,prime);   %y^2*z
quad21 = FF(P2.value(3,:),prime)*FF(lin_vars,prime);   %y*z^2
cube_2 = FF(P2.value(4,:),prime)*FF(lin_vars,prime);  %z^3

identities = [cube_1*cube_2 - quad12*quad21;...
    (cube_1)*FF(lin_vars(5),prime) - (quad21)*FF(lin_vars(4),prime);...                % (z^3)*y = (z^2*y)*z
    (cube_1)*FF(lin_vars(5)^2,prime) - (quad12)*FF(lin_vars(4),prime);...
    (cube_2)*FF(lin_vars(4),prime) - (quad21)*FF(lin_vars(5),prime);...                % (y^2*z)*z = (z^2*y)*y
    (cube_2)*FF(lin_vars(4)^2,prime) - (quad12)*FF(lin_vars(5)^2,prime);...                     % (y^2*z) = (y*z)*y
    (quad12)*FF(lin_vars(5),prime) - (quad21)*FF(lin_vars(4),prime)];% (y*z)*z = (z^2*y)
for i=1:6
    identities(i) = collect(identities(i));
end
old_vars = [lin_vars(4)^3,lin_vars(5)^3,lin_vars(4)^2*lin_vars(5),lin_vars(4)*lin_vars(5)^2];
new_vars = [cube_1.value,cube_2.value,quad12.value,quad21.value];
for i=1:6
    identities(i) = expand(subs(expand(identities(i)),old_vars,new_vars));
end
u = sym("u","real");
v = sym("v","real");
w = sym("w","real");
for i=1:6
    identities(i) = subs(identities(i),all_lin_vars(1:3),[u;v;w]);
end
eq = [identities(1).value,identities(2).value,identities(3).value,identities(4).value,identities(5).value,identities(6).value];
eq = simplify(eq);
[A_,b] = equationsToMatrix(eq,[u;v;w;lin_vars(4:5)]);
A2 = FF([A_,-b],prime);
coef = coeffs(det(A2),p_var,"All");
p_root = get_gf_root(coef,prime);
result = [];
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars(4)),find(vars==lin_vars(5))]);
for i=1:numel(p_root)
    if p_root(i) == 0
        p_root(i) = p_root(i);
    end
    M = subs(A2,p_var,p_root(i));
    rank(M)
    % cur_result = solve_subsystem4C3(M,p_root(i),idx,plot_subspace=0);
    rM = rref(M);
    col = zeros(1,rank(rM));
    for j=1:rank(rM)
        col(j) = find(rM(j,:),1,"first");
    end
    m = size(rM,1);
    vec = 1:m;
    r = rM(1:rank(M),setdiff(vec,col));
    cur_result = [];
    disp(join(string(col),""))
    switch join(string(col),"")
        case "1"
            [u_root,v_root] = rank1_4C3_1_fp(r,prime);
        case "2"
            [u_root,v_root] = rank1_4C3_2_fp(r,prime);
        case "3"
            [u_root,v_root] = rank1_4C3_3_fp(r,prime);
        case "4"
            [u_root,v_root] = rank1_4C3_4_fp(r,prime);
        case "12"
            [u_root,v_root] = rank2_4C3_12_fp(r,prime);
        case "13"
            [u_root,v_root] = rank2_4C3_13_fp(r,prime);
        case "14"
            [u_root,v_root] = rank2_4C3_14_fp(r,prime);
        case "23"
            [u_root,v_root] = rank2_4C3_23_fp(r,prime);
        case "24"
            [u_root,v_root] = rank2_4C3_24_fp(r,prime);
        case "34"
            [u_root,v_root] = rank2_4C3_34_fp(r,prime);
        case "123"
            [u_root,v_root] = rank3_4C3_123_fp(r,prime);
        case "124"
            [u_root,v_root] = rank3_4C3_124_fp(r,prime);
        case "134"
            [u_root,v_root] = rank3_4C3_134_fp(r,prime);
        case "234"
            [u_root,v_root] = rank3_4C3_234_fp(r,prime);
        case "1234"
            [u_root,v_root] = rank4_4C3_1234_fp(r,prime);
        case "1245"
            [u_root,v_root] = rank4_4C3_1234_fp(r,prime);
        case "1345"
            [u_root,v_root] = rank4_4C3_1234_fp(r,prime);
        case "2345"
            [u_root,v_root] = rank4_4C3_1234_fp(r,prime);
        case "12345"
            [u_root,v_root] = rank5_4C3_12345_fp(r,prime);
        otherwise
            continue
            % fprintf("System of equations is singular.\n");
    end
    if ~isempty(u_root)
        cur_result=[repmat(p_root,numel(u_root),1),reshape(u_root,[],1),reshape(v_root,[],1)];
    end
    if ~isempty(cur_result)
        result = [result;cur_result];
    end
end
if ~isempty(result)
    result = result(:,idx);
end
completion_time = toc(t1);
fprintf("Algorithm completed in %.2fs.\n",completion_time);
equations = c*var_vec;
if numel(result) > 0
    print_solutions(result,equations,x,y,z)
end
end