function [result] = E3C3(C,opt)
%E3C3 Main function for the 3C3 implementation.
arguments
    C (3,20) {mustBeReal};
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 1;
    opt.plot_surf {mustBeInRange(opt.plot_surf,0,1)} = 1;
    opt.plot_subspace {mustBeInRange(opt.plot_subspace,0,1)} = 0;
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
    opt.tolerance (1,1) {mustBeReal,mustBePositive} = 10^-10;
    opt.error (1,1) {mustBeNumericOrLogical} = true;
    opt.log_db (1,1) {mustBeNumericOrLogical} = false;
    opt.complex (1,1) {mustBeNumericOrLogical} = false;
end
opt = check_toolboxes(opt);
t1 = tic;
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
vars = [x,y,z];
var_vec = [x.^3,x.^2.*y,x.^2.*z,x.*y.^2,...
    x.*y.*z,x.*z.^2,y.^3,y.^2.*z,y.*z.^2,z.^3,x.^2,x.*y,x.*z,y.^2,y.*z,z.^2,x,y,z,1]';
result = [];
% Raise error or output warning if split_matrices fails
try
    [lin_vars,p_var,P2,prc] = split_matrices_3C3(C,x,y,z,verbose=opt.verbose,error=opt.error);
catch exception
    if opt.error
        rethrow(exception)
    else
        warning(exception.identifier,"%s",exception.message)
        return
    end
end

A = substitute_identities_3C3(P2,lin_vars,prc);
% Compute p_var solutions
coef = coeffs(numden(det(A)),p_var,"All");
p_root = roots(coef);
if ~opt.complex
    p_root = p_root(abs(imag(p_root))<10^-10);
else
    p_root(abs(imag(p_root))<10^-10) = p_root(abs(imag(p_root))<10^-10);
end
p_root = unique(p_root);
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars(5)),find(vars==lin_vars(6))]);
for i=1:numel(p_root)
    M = subs(A,p_var,p_root(i));
    if opt.verbose > 1
        disp(M)
        disp(rref(double(M),opt.tolerance))
    end
    cur_result = solve_subsystem_3C3(M,p_root(i),plot_subspace=opt.plot_subspace,tolerance=opt.tolerance,verbose=opt.verbose,complex=opt.complex);
    if ~isempty(cur_result)
        result = [result;cur_result]; %#ok<AGROW>
    end
end
if ~isempty(result)
    result = result(:,idx);
end
completion_time = toc(t1);
fprintf("Algorithm completed in %.2fs.\n",completion_time);
if opt.plot_surf
    plot_and_color_implicit(result,3,C,tolerance=opt.tolerance)
end
equations = C*var_vec;
% Output solutions to console
if ~isempty(result) && opt.verbose
    print_solutions(result,equations,x,y,z)
end
% Log run details to database
if opt.log_db
    log_to_db(C,result,completion_time,opt.tolerance,0)
end
end