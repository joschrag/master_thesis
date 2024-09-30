function [result] = E3C3_Fp(C,prime,opt)
%E3C3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    C (3,20) {mustBeInteger};
    prime (1,1) {mustBePrime};
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 1;
    opt.error (1,1) {mustBeNumericOrLogical} = true;
    opt.log_db (1,1) {mustBeNumericOrLogical} = false;
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
    [lin_vars, p_var, P2] = split_matrices_3C3_Fp(C,prime,x,y,z,verbose=opt.verbose);
catch exception
    if opt.error
        rethrow(exception)
    else
        warning(exception.identifier,"%s",exception.message)
        return
    end
end

A = substitute_identities_3C3_Fp(P2,lin_vars,prime);
% Compute p_var solutions
pol = coeffs(det(A),p_var,"All");
% Since we cant use rank() for polynomial matrices over Fp, we check wether the
% determinant is the zero polynomial
if isempty(pol)
    if opt.error
        error("Polynomial matrix has non-full rank!")
    else
        warning("Polynomial matrix has non-full rank!")
        return
    end
end
p_root = get_gf_root(pol,prime);
equations = C*var_vec;
% Determine variable indices to reorder solutions
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars(5)),find(vars==lin_vars(6))]);
for i=1:numel(p_root)
    M = FF(subs(A,p_var,p_root(i)),prime);
    cur_result = solve_subsystem_3C3_Fp(M,p_root(i),prime,verbose=opt.verbose);
    if ~isempty(cur_result)
        result = [result;cur_result]; %#ok<AGROW>
    end
end
if ~isempty(result)
    result = uint64(result(:,idx));
end
completion_time = toc(t1);
fprintf("Algorithm completed in %.2fs.\n",completion_time);
% Output solutions to console
if ~isempty(result) && opt.verbose
    print_solutions(result,equations,x,y,z,prime)
end
% Log run details to database
if opt.log_db
    log_to_db(C,result,completion_time,0,prime);
end
end