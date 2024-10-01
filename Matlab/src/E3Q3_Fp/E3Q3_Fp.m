function [result] = E3Q3_Fp(C,prime,opt)
%E3Q3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    C (:,10) {mustBeInteger};
    prime (1,1) {mustBePrimeOrZero,mustBePositive};
    opt.verbose (1,1) {mustBeInteger, mustBeInRange(opt.verbose,0,2)} = 1;
    opt.error (1,1) {mustBeNumericOrLogical} = true;
    opt.log_db (1,1) {mustBeNumericOrLogical} = false;
end
opt = check_toolboxes(opt);
t1 = tic;
x = sym("x","integer");
y = sym("y","integer");
z = sym("z","integer");
vars = [x,y,z];
v = [x^2,x*y,x*z,y^2,y*z,z^2,x,y,z,1]';
result = [];
% Raise error or output warning if split_matrices fails
try
    [lin_vars, p_var, P2] = split_matrices_E3Q3_Fp(C,prime,x,y,z,verbose=opt.verbose);
catch exception
    if opt.error
        rethrow(exception)
    else
        warning(exception.identifier,"%s",exception.message)
        return
    end
end
A = substitute_identities_E3Q3_Fp(P2,lin_vars,prime);
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
if opt.verbose > 1
    fprintf("p_root: %s\n",mat2str(p_root))
end
% Determine variable indices to reorder solutions
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars(1)),find(vars==lin_vars(2))]);
result = [];
equations = C*v;
for i = 1:numel(p_root)
    p_0 = p_root(i);
    M = subs(A,p_var,p_0);
    cur_result = solve_subsystem_E3Q3_Fp(M,p_0,prime,equations,p_var,lin_vars(1:2),verbose=opt.verbose);
    if ~isempty(cur_result)
        result = [result;cur_result]; %#ok<AGROW>
    end
end
if ~isempty(result)
    result = uint64(result(:,idx));
end
completion_time = toc(t1);
if opt.verbose > 0
    fprintf("Algorithm completed in %.2fs.\n",completion_time);
end
% Output solutions to console
if ~isempty(result) && opt.verbose
    print_solutions(result,equations,x,y,z,prime)
end
% Log run details to database
if opt.log_db
    log_to_db(C,result,completion_time,0,prime);
end
end