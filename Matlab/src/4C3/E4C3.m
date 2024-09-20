function [result] = E4C3(c,opt)
arguments       %1  2   3   4   5   6   7   8   9  10  11  12   13 14 15  16   17  18  19  20
    c (4,20) = [0,-3,-2,0,0,0,-1,0,0,0,0,6,2,0,0,0,-15,-5,2,-5;...
        0,2,-1,0,0,0,0,-1,0,0,0,2,4,0,0,0,10,-8,2,-2;...
        -1,2,1,0,0,0,0,0,-1,0,0,2,-8,0,0,0,3,7,13,-4;...
        1,2,-2,0,0,0,0,0,0,-1,0,-2,1,0,0,0,-7,-4,-5,-5;...
        ];
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 0;
    opt.plot_surf {mustBeInRange(opt.plot_surf,0,1)} = 1;
    opt.plot_subspace {mustBeInRange(opt.plot_subspace,0,1)} = 1;
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
    opt.tolerance (1,1) {mustBeReal,mustBePositive} = 10^-10;
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
conds = zeros(1,3);
for i=1:3
    Q{i} = -[c(:,q_idx(i,1)),c(:,q_idx(i,2)),c(:,q_idx(i,3)),c(:,q_idx(i,4))];
    ranks(i) = rank(Q{i});
    conds(i) = real(cond(Q{i}));
end
if any(ranks==4)
    [M,I] = min(conds,[],"omitnan");
else
    M = Inf;
end

if ~isinf(M)
    var_idx = setdiff(1:3,I);
all_lin_vars = [vars(var_idx(1)).^2,vars(var_idx(1))*vars(var_idx(2)),vars(var_idx(2)).^2,vars(setdiff(1:3,I)),1]';
    if opt.verbose > 0
        fprintf("Using P(%s)\n",string(vars(I)))
    end
    p_var = vars(I);
    p_var_pow = [1;p_var;p_var^2;p_var^3];
    P = sym.zeros(4,6);
    for j=1:6
        P(:,j) = c(:,complete_idx{I,j})*p_var_pow(numel(complete_idx{I,j}):-1:1);
    end
    lin_vars = all_lin_vars;
    P2 = Q{I}\P;
else
    result = [];
    return
end

if ~all(all(P2(:,1:1:3)==sym(0)))
    % warning("Substitution not possible!")
    result = [];
    return
end

cube_1 = P2(1,:)*lin_vars;   %y^3
quad12 = P2(2,:)*lin_vars;   %y^2*z
quad21 = P2(3,:)*lin_vars;   %y*z^2
cube_2 = P2(4,:)*lin_vars;  %z^3

identities = [cube_1*cube_2 - quad12*quad21;...
    (cube_1)*lin_vars(5) - (quad21)*lin_vars(4);...                % (z^3)*y = (z^2*y)*z
    (cube_1)*lin_vars(5)^2 - (quad12)*lin_vars(4);...
    (cube_2)*lin_vars(4) - (quad21)*lin_vars(5);...                % (y^2*z)*z = (z^2*y)*y
    (cube_2)*lin_vars(4)^2 - (quad12)*lin_vars(5)^2;...                     % (y^2*z) = (y*z)*y
    (quad12)*lin_vars(5) - (quad21)*lin_vars(4)];                      % (y*z)*z = (z^2*y)
sub_new = collect(identities);
old_vars = [lin_vars(4)^3,lin_vars(5)^3,lin_vars(4)^2*lin_vars(5),lin_vars(4)*lin_vars(5)^2];
new_vars = collect([cube_1,cube_2,quad12,quad21],lin_vars(4:5));
sub_new = expand(subs(expand(sub_new),old_vars,new_vars));
u = sym("u","real");
v = sym("v","real");
w = sym("w","real");
sub_final = subs(sub_new,all_lin_vars(1:3),[u;v;w]);
[A_,b] = equationsToMatrix(sub_final,[u;v;w;lin_vars(4:5)]);
A2 = [A_,-b];
coef = coeffs(det(A2),p_var,"All");
p_root = roots(coef);
p_root = unique(p_root(abs(imag(p_root))<10^-10));
result = [];
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars(4)),find(vars==lin_vars(5))]);
for i=1:numel(p_root)
    if p_root(i) == 0
        p_root(i) = p_root(i);
    end
    M = subs(A2,p_var,p_root(i));
    cur_result = solve_subsystem4C3(M,p_root(i),idx,plot_subspace=0);
    if ~isempty(cur_result)
        result = [result;cur_result];
    end
end
completion_time = toc(t1);
fprintf("Algorithm completed in %.2fs.\n",completion_time);
if ~isempty(result) && opt.plot_surf
    plot_and_color_implicit(result,m,c)
end
equations = c*var_vec;
if numel(result) > 0
    print_solutions(result,equations,x,y,z)
end
end