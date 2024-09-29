function [result] = E4C3_Fp(C,prime,opt)
arguments       %1  2   3   4   5   6   7   8   9  10  11  12   13 14 15  16   17  18  19  20
    C (4,20) = [0,-3,-2,0,0,0,-1,0,0,0,0,6,2,0,0,0,-15,-5,2,-5;...
        0,2,-1,0,0,0,0,-1,0,0,0,2,4,0,0,0,10,-8,2,-2;...
        -1,2,1,0,0,0,0,0,-1,0,0,2,-8,0,0,0,3,7,13,-4;...
        1,2,-2,0,0,0,0,0,0,-1,0,-2,1,0,0,0,-7,-4,-5,-5;...
        ];
    prime (1,1) {mustBePrime} = nextprime(6);
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 1;
    opt.error (1,1) {mustBeNumericOrLogical} = true;
    opt.log_db (1,1) {mustBeNumericOrLogical} = true;
end
result=[];
t1 = tic;
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
vars = [x,y,z];
var_vec = [x.^3,x.^2.*y,x.^2.*z,x.*y.^2,...
    x.*y.*z,x.*z.^2,y.^3,y.^2.*z,y.*z.^2,z.^3,x.^2,x.*y,x.*z,y.^2,y.*z,z.^2,x,y,z,1]';
try
    [lin_vars, p_var, P2] = split_matrices_4C3_Fp(C,prime,x,y,z,verbose=opt.verbose);
catch exception
    if opt.error
        rethrow(exception)
    else
        warning(exception.identifier,"%s",exception.message)
        return
    end
end

A = substitute_identities_4C3_Fp(P2,lin_vars,prime);
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
% Determine variable indices to reorder solutions
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars(4)),find(vars==lin_vars(5))]);
for i=1:numel(p_root)
    M = subs(A,p_var,p_root(i));
    cur_result = solve_subsystem_4C3_Fp(M,p_root(i),prime);
    if ~isempty(cur_result)
        result = [result;cur_result];
    end
end
if ~isempty(result)
    result = result(:,idx);
end
completion_time = toc(t1);
fprintf("Algorithm completed in %.2fs.\n",completion_time);
equations = C*var_vec;
% Output solutions to console
if ~isempty(result)
    print_solutions(result,equations,x,y,z,prime)
end
% Log run details to database
% for i=0:prime-1
%     for j=0:prime-1
%         for k=0:prime-1
%             if mod(subs(equations,[x,y,z],[i,j,k]),prime)==0
%                 disp([i,j,k])
%             end
%         end
%     end
% end
end