function [result] = E3C3_Fp(c,prime,opt)
%E3C3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    c (3,20) = [     3,     0,     0 ,    0,     0,     0,    -1,     0 ,    0 ,    0,     4,     0,     0,  0,     0,     0,     2,     0,     0,     2;...
     2,     0 ,    0 ,    0,     0,     0,     0,    -1,     0,     0,     3,     0,     0,  0,     0,     0,     4,     0,     0,     2;...
     1,     0,     0,     0 ,    0 ,    0,     0,     0 ,    0,     0,     0,     0,     0,  0,    -1,     0,     3,     0,     0,     4];
    prime (1,1) {mustBeInteger,mustBePositive} = nextprime(4);
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
complete_idx = {9,10,[4,14],[6,16],[2,12,18],[3,13,19],[1,11,17,20];...
    6,10,[2,11],[9,16],[4,12,17],[8,15,19],[7,14,18,20];...
    4,7,[3,11],[8,14],[6,13,17],[9,15,18],[10,16,19,20]};
q_idx = [7,8,5,15;1,3,5,13;1,2,5,12];    %indices of cubic monomials of lin_vars
assert(all(mod(c(:,5),prime)==zeros(3,1)))
Q = cell(1,3);
ranks = zeros(1,3);
for i=1:3
    Q{i} = FF(-[c(:,q_idx(i,1)),c(:,q_idx(i,2)),c(:,q_idx(i,3))*x+c(:,q_idx(i,4))],prime);
    ranks(i) = rank(Q{i});
end
valid_P2 = false;

if any(ranks==3)
    I = find(ranks==3);
    for idx = I
        var_idx = setdiff(1:3,idx);
        all_lin_vars = [vars(var_idx(1))*vars(var_idx(2))^2,vars(var_idx(2)).^3,vars(var_idx(1))^2,vars(var_idx(2)).^2,vars(setdiff(1:3,I)),1]';
        if opt.verbose > 0
            fprintf("Using P(%s)\n",string(vars(idx)))
        end
        p_var = vars(idx);
        p_var_pow = [1;p_var;p_var^2;p_var^3];
        P = sym.zeros(3,7);
        for j=1:7
            P(:,j) = c(:,complete_idx{idx,j})*p_var_pow(numel(complete_idx{idx,j}):-1:1);
        end
        lin_vars = all_lin_vars;
        P2 = Q{idx}^(-1)*FF(P,prime);
        if all(all(P2.value(:,1:6)==sym(0)))
            valid_P2 = true;
            break
        end
    end
else
    warning("All matrices singular!")
    return
end

if ~valid_P2
    warning("Substitution not possible!")
    result = [];
    return
end

cube = FF(P2.value(1,:)*lin_vars,prime);  %y^3
quad = FF(P2.value(2,:)*lin_vars,prime);   %y^3
mixed = FF(P2.value(3,:)*lin_vars,prime);   %y^3

identities = [...
    cube*FF(lin_vars(6),prime)-quad*FF(lin_vars(5),prime);...
    cube*FF(lin_vars(6),prime) - mixed*FF(lin_vars(5)^2,prime);...
    quad*FF(lin_vars(6),prime)-mixed^2;...
    cube*FF(lin_vars(6)^3,prime)-mixed^3;...
    quad-mixed*FF(lin_vars(5),prime);...
    % cube*FF(lin_vars(6)^3,prime)-quad*FF(lin_vars(5)*lin_vars(6)^2,prime);...
    cube*FF(lin_vars(6)^2,prime)-quad*mixed;...
    ];
for i=1:6
    identities(i) = collect(identities(i));
end
old_vars = [lin_vars(5)^3,lin_vars(5)^2*lin_vars(6),lin_vars(5)*lin_vars(6)];
new_vars = collect([cube.value,quad.value,mixed.value],lin_vars(5:6));
for i=1:6
    identities(i) = expand(subs(expand(identities(i)),old_vars,new_vars));
end
u = sym("u","real");
v = sym("v","real");
w = sym("w","real");
for i=1:6
    identities(i) = subs(identities(i),all_lin_vars(2:4),[u;v;w]);
end
eq = [identities(1).value,identities(2).value,identities(3).value,identities(4).value,identities(5).value,identities(6).value];
[A_,b] = equationsToMatrix(eq,[u;v;w;lin_vars(5:6)]);
A2 = [A_,-b];
coef = coeffs(det(A2),p_var,"All");
p_root = get_gf_root(coef,prime)
result = [];
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars(5)),find(vars==lin_vars(6))]);
for i=1:numel(p_root)
    M = FF(subs(A2,p_var,p_root(i)),prime);
    cur_result = solve_subsystem3C3_Fp(M,p_root(i),idx,prime);
    if ~isempty(cur_result)
        result = [result;cur_result];
    end
end
completion_time = toc(t1);
fprintf("Algorithm completed in %.2fs.\n",completion_time);
equations = c*var_vec;
if numel(result) > 0
    print_solutions(result,equations,x,y,z,prime)
end
% disp(equations)
for i=0:prime-1
    for j=0:prime-1
        for k=0:prime-1
            if mod(subs(c*var_vec,[x,y,z],[i,j,k]),prime)==0
                disp([i,j,k])
            end
        end
    end
end
end



