function [result] = E3Q3_Fp(c,prime,opt)
%E3Q3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    c (:,10) {mustBeInteger} = [1, 2, 1, 0, 1, 1, 2, 0, 0, -4;...
        1, 1, 3, 2, 1, -1, 2, 0, 2, -5;...
        1, -1, 1, 1, 0, 1, -1, 0, 3, -6;
        ];
    prime (1,1) {mustBeInteger} = nextprime(6);
    opt.verbose (1,1) {mustBeInteger, mustBeInRange(opt.verbose,0,2)} = 1;
    opt.error (1,1) {mustBeNumericOrLogical} = true;
    opt.log_db (1,1) {mustBeNumericOrLogical} = true;
end
t1 = tic;
x = sym("x","integer");
y = sym("y","integer");
z = sym("z","integer");
vars = [x,y,z];
v = [x^2,x*y,x*z,y^2,y*z,z^2,x,y,z,1]';
[lin_vars, p_var, P2] = split_matrices_Fp(c,prime,size(c,1),x,y,z,verbose=opt.verbose,error=opt.error);

quad_1 = FF(P2.value(1,:),prime)*lin_vars;
quad_2 = FF(P2.value(3,:),prime)*lin_vars;
mixed = FF(P2.value(2,:),prime)*lin_vars;

identities = [(quad_1)*FF(lin_vars.value(2),prime) - (mixed)*FF(lin_vars.value(1),prime);...
    (mixed)*FF(lin_vars.value(2),prime) - (quad_2)*FF(lin_vars.value(1),prime);...
    (mixed)*(mixed) - quad_1*quad_2];
sub1 = FF.empty(3,0);
sub2 = FF.empty(3,0);
sub3 = FF.empty(3,0);

for i = 1:3
    sub1(i) = collect(identities(i));
end
old_vars = [lin_vars.value(1)^2,...
    lin_vars.value(2)^2,...
    lin_vars.value(1)*lin_vars.value(2)];
new_vars = [quad_1.value,quad_2.value,mixed.value];
for i = 1:3
    sub2(i) = subs(sub1(i),old_vars,new_vars);
end
for i = 1:3
    ex = expand(sub2(i));
    sub3(i) = subs(ex,old_vars,new_vars);
end
eq = [sub3(1).value,sub3(2).value,sub3(3).value];
eq = simplify(eq);
[A_,b] = equationsToMatrix(eq,lin_vars.value(1:2));
A2 = FF([A_,-b],prime);
pol = coeffs(det(A2),p_var,"All");
p_root = get_gf_root(pol,prime);
if numel(p_root) == 0
    fprintf("No solutions to the equations.\n")
    result = [];
    return;
end
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars.value(1)),find(vars==lin_vars.value(2))]);
result = [];
for i = 1:numel(p_root)
    p_0 = p_root(i);
    M = subs(A2,p_var,p_0);
    if M.value == zeros(3)
        fprintf("Matrix is singular.\n");
        continue;
    end
    rM = rref(M);
    m = zeros(1,rank(rM));
    for j=1:rank(rM)
        m(j) = find(rM(j,:),1,"first");
    end
    r = vpa(rM(1:rank(rM),setdiff(1:3,m)));
    switch join(string(m),"")
        case "1"
            r_root = (0:prime-1)';
            q_root = FF(-r(1).*r_root-r(2),prime).value;
        case "2"
            q_root = (0:prime-1)';
            r_root = repmat(FF(-r(2),prime).value,size(q_root));
        case "12"
            o_sols = FF(-r',prime).value;
            q_root = o_sols(1);
            r_root = o_sols(2);
        otherwise
            fprintf("Matrix equations are inconsistent!");
            continue;
    end
    if ~isempty(q_root)
        for root=[q_root,r_root]'
            result=[result;[p_0,root']];
        end
    end
end
result = uint64(result(:,idx));
completion_time = toc(t1);
if opt.log_db
    log_to_db(c,result,completion_time,0,prime);
end
equations = c*v;
if numel(result) > 0
    print_solutions(result,equations,x,y,z,prime)
end
end

