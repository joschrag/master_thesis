function [result] = E5C3(c,opt)
%E3Q3_CUBIC Summary of this function goes here
%   Detailed explanation goes here
arguments       %1  2   3   4   5   6   7   8   9  10  11  12   13 14 15  16   17  18  19  20
    c (5,20) = [1,-3,-2,0,0,0,-1,4,-1,-2,2,6,2,-1,-4,5,-1,-5,2,-2;...
        1,2,-1,-3,0,-1,0,1,-2,2,-4,2,4,3,4,-5,2,-8,2,1;...
        -2,2,1,-4,0,2,-5,0,2,2,1,2,-8,22,-6,-8,9,-20,13,-3;...
        -2,2,-2,-1,0,2,2,-2,-2,2,6,-2,1,-4,9,-5,-4,-4,-5,9;...
        0,-2,0,-1,0,0,-3,0,4,0,-2,6,1,6,-8,-2,3,2,3,-7;...
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
[lin_vars,p_var,P2] = split_matrices(c,m,x,y,z,verbose=opt.verbose,error=opt.error);

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
        if opt.error
            error("Precision of root is too low!")
        else
            warning("Precision of root is too low!")
        end
    end
    cur_result = solve_subsystem5C3(M,p_root(i),idx,plot_subspace=opt.plot_subspace);
    if ~isempty(cur_result)
        result = [result;cur_result];
    end
end
completion_time = toc(t1);
if opt.plot_surf
    plot_and_color_implicit(result, m, c);
end
if opt.log_db
    log_to_db(c,result,completion_time,opt.tolerance,0)
end
equations = c*var_vec;
if numel(result) > 0
    print_solutions(result,equations,x,y,z)
end
end

