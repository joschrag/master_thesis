function [result] = E4C3(c,opt)
arguments       %1  2   3   4   5   6   7   8   9  10  11  12   13 14 15  16   17  18  19  20
    c (4,20) = [0,-3,-2,0,0,0,-1,0,0,0,0,6,2,0,0,0,-15,-5,2,-5;...
        0,2,-1,0,0,0,0,-1,0,0,0,2,4,0,0,0,10,-8,2,-2;...
        -1,2,1,0,0,0,0,0,-1,0,0,2,-8,0,0,0,3,7,13,-4;...
        1,2,-2,0,0,0,0,0,0,-1,0,-2,1,0,0,0,-7,-4,-5,-5;...
        ];
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 1;
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
try
    [lin_vars,p_var,P2] = split_matrices_4C3(c,x,y,z,verbose=opt.verbose);
catch exception
    if opt.error
        rethrow(exception)
    else
        warning(exception.identifier,"%s",exception.message)
        return
    end
end

A = substitute_identities_4C3(P2,lin_vars);
coef = coeffs(det(A),p_var,"All");
p_root = roots(coef);
p_root = unique(p_root(abs(imag(p_root))<10^-10));
result = [];
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars(4)),find(vars==lin_vars(5))]);
for i=1:numel(p_root)
    if p_root(i) == 0
        p_root(i) = p_root(i);
    end
    M = subs(A,p_var,p_root(i));
    cur_result = solve_subsystem_4C3(M,p_root(i),idx,plot_subspace=0);
    if ~isempty(cur_result)
        result = [result;cur_result];
    end
end
completion_time = toc(t1);
fprintf("Algorithm completed in %.2fs.\n",completion_time);
if  opt.plot_surf % ~isempty(result) && opt.plot_surf
    plot_and_color_implicit(result,4,c)
end
equations = c*var_vec;
if numel(result) > 0
    print_solutions(result,equations,x,y,z)
end
end