function [result] = E3Q3(C,opt)
%E3Q3 Implements the E3Q3 algorithm over the real numbers.
arguments
    C (3,10) {mustBeReal} = [...
        1,1,1,2,0,1,1,1,1,-10;...
        1,0,1,1,-1,3,7,0,0,-10;...
        1,0,0,1,1,1,-15,0,0,-10;...
        ];
    opt.verbose (1,1) {mustBeInteger, mustBeInRange(opt.verbose,0,2)} = 1;
    opt.plot (1,1) {mustBeNumericOrLogical} = true;
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
    opt.tolerance (1,1) {mustBeReal,mustBePositive} = 10^-10;
    opt.error (1,1) {mustBeNumericOrLogical} = true;
    opt.log_db (1,1) {mustBeNumericOrLogical} = false;
    opt.show_lines (1,1) {mustBeNumericOrLogical} = false;
end
opt = check_toolboxes(opt);
t1 = tic;
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
vars = [x,y,z];
v = [x^2,x*y,x*z,y^2,y*z,z^2,x,y,z,1]';
result = [];
% Raise error or output warning if split_matrices fails
try
    [lin_vars, p_var, P2] = split_matrices_E3Q3(C, x, y, z,verbose=opt.verbose);
catch exception
    if opt.error
        rethrow(exception)
    else
        warning(exception.identifier,"%s",exception.message)
        return
    end
end
A = substitute_identities_E3Q3(P2,lin_vars);
% Compute p_var solutions
coef = coeffs(det(A),p_var,"All");
p_root = unique(roots(coef));
p_root = p_root(imag(p_root)==0);
% Determine variable indices to reorder solutions
[~,idx] = sort([find(vars==p_var),find(vars==lin_vars(1)),find(vars==lin_vars(2))]);
result = [];
equations = C*v;
for i =1:numel(p_root)
    M = subs(A,p_var,p_root(i));
    if rank(double(M),opt.tolerance) == 3
        if opt.error
            error("Precision of root is too low!")
        else
            warning("Precision of root is too low!")
            continue
        end
    end
    if opt.verbose > 1
        disp(M)
        disp(rref(double(M),opt.tolerance))
    end
    cur_result = solve_subsystem_E3Q3(M,p_root(i),equations,p_var,lin_vars,...
        tolerance=opt.tolerance,show_lines=opt.show_lines);
    if ~isempty(cur_result)
        result = [result;cur_result]; %#ok<AGROW>
    end
end
if ~isempty(result)
    result = result(:,idx);
end
completion_time = toc(t1);
% Plot solutions and surfaces
if opt.plot
    ax = gca();
    s = scatter3(ax,result(:,1),result(:,2),result(:,3),30,"k","filled");
    set(s,"DisplayName","solution")
    hold(ax,"on")
    min_axes = min(result);
    max_axes = max(result);
    min_dist = 0.1;
    axis_vec =  zeros(1,6);
    for i=1:3
        extra = abs(max(max_axes(i)-min_axes(i),min_dist))*0.35;
        axis_vec(2*i-1:2*i) = [min_axes(i)-extra,max_axes(i)+extra];
    end
    plot_and_color(C,verbose=opt.verbose,plotRange=opt.plotRange,tolerance=opt.tolerance)
    axis(double(axis_vec))
end
fprintf("Algorithm completed in %.2fs.\n",completion_time);
% Output solutions to console
if ~isempty(result) && opt.verbose
    print_solutions(result,equations,x,y,z)
end
% Log run details to database
if opt.log_db
    log_to_db(C,result,completion_time,opt.tolerance,0)
end
end