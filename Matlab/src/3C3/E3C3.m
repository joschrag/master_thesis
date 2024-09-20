function [result] = E3C3(c,opt)
%E3C3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    c (3,20) = [5,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,-15,0,0,-5;...
                0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,10,0,0,-2;...
                -1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,3,0,0,-4;...
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
complete_idx = {9,10,[4,14],[6,16],[2,12,18],[3,13,19],[1,11,17,20];...
    6,10,[2,11],[9,16],[4,12,17],[8,15,19],[7,14,18,20];...
    4,7,[3,11],[8,14],[6,13,17],[9,15,18],[10,16,19,20]};
q_idx = [7,8,5,15;1,3,5,13;1,2,5,12];    %indices of cubic monomials of lin_vars
assert(all(c(:,5)==zeros(3,1)))
Q = cell(1,3);
ranks = zeros(1,3);
conds = zeros(1,3);
for i=1:3
    Q{i} = -[c(:,q_idx(i,1)),c(:,q_idx(i,2)),c(:,q_idx(i,3))*x+c(:,q_idx(i,4))];
    ranks(i) = rank(Q{i});
    conds(i) = real(cond(Q{i}));
end
if any(ranks==3)
    [M,I] = min(conds,[],"omitnan");
else
    M = Inf;
end
if ~isinf(M)
    var_idx = setdiff(1:3,I);
    all_lin_vars = [vars(var_idx(1))*vars(var_idx(2)).^2,vars(var_idx(2)).^3,vars(var_idx(1))^2,vars(var_idx(2)).^2,vars(setdiff(1:3,I)),1]';
    if opt.verbose > 0
        fprintf("Using P(%s)\n",string(vars(I)))
    end
    p_var = vars(I);
    p_var_pow = [1;p_var;p_var^2;p_var^3];
    P = sym.zeros(3,7);
    for j=1:7
        P(:,j) = c(:,complete_idx{I,j})*p_var_pow(numel(complete_idx{I,j}):-1:1);
    end
    lin_vars = all_lin_vars;
    P2 = Q{I}\P;
else
    warning("All Matrices are singular!")
    result = [];
    return
end
if ~all(all(P2(:,1:6)==sym(0)))
    warning("Substitution not possible!")
    result = [];
    return
end
cube = P2(1,:)*lin_vars;  %y^3
quad = P2(2,:)*lin_vars;   %y^3
mixed = P2(3,:)*lin_vars;   %y^3

identities = [...
    cube*lin_vars(6)-quad*lin_vars(5);...
    cube*lin_vars(6) - mixed*lin_vars(5)^2;...
    quad*lin_vars(6)-mixed^2;...
    cube*lin_vars(6)^3-mixed^3;...
    quad-mixed*lin_vars(5);...
    % cube*lin_vars(6)^3-quad*lin_vars(5)*lin_vars(6)^2;...
    cube*lin_vars(6)^2-quad*mixed;...
    ];
sub_new = collect(identities);
old_vars = [lin_vars(5)^3,lin_vars(5)^2*lin_vars(6),lin_vars(5)*lin_vars(6)];
new_vars = collect([cube,quad,mixed],lin_vars(5:6));
sub_new = expand(subs(expand(sub_new),old_vars,new_vars));
u = sym("u","real");
v = sym("v","real");
w = sym("w","real");
sub_final = subs(sub_new,lin_vars(2:4),[u;v;w]);
[A_,b] = equationsToMatrix(sub_final,[u;v;w;lin_vars(5:6)]);
A2 = [A_,-b];
coef = coeffs(det(A2),p_var,"All");
p_root = roots(coef);
p_root = unique(p_root(abs(imag(p_root))<10^-10));
result = [];
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars(5)),find(vars==lin_vars(6))]);
for i=1:numel(p_root)
    M = subs(A2,p_var,p_root(i));
    cur_result = solve_subsystem3C3(M,p_root(i),idx,plot_subspace=opt.plot_subspace,tolerance=opt.tolerance);
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

