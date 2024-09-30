function [result] = solve_subsystem_E3Q3(M,p_root,equations,p_var,lin_vars,opt)
%SOLVE_SUBSYSTEM_E3Q3 Solve the resulting system of equations for all cases.
%   For all cases of rref(M) this function solves the system of equations and
%   returns the solution.
arguments
    M (3,3) sym {mustBeReal};
    p_root (1,1) sym {mustBeReal};
    equations (3,1) sym {mustBeReal};
    p_var (1,1) sym;
    lin_vars (3,1) sym;
    opt.show_lines (1,1) {mustBeNumericOrLogical} = false;
    opt.tolerance (1,1) {mustBeReal,mustBePositive} = 10^-10;
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 1;
end
result = [];
if rank(M)==0
    warning("Zero matrix for p_root %i",p_root)
    return
end
if rank(double(M),opt.tolerance) == 3
    warning("Precision of root %f is too low!",p_root)
    return
end
switch rank(double(M),opt.tolerance)
    % If rank(M)=1 go through different cases and solve based on rref matrix
    case 1
        rM = rref(double(M), opt.tolerance);
        if rank(M) ~= rank(rM)
            warning("Numerical error calculating rref matrix for p_root %i",p_root)
            return
        end
        col = zeros(1,rank(rM));
        for j=1:rank(rM)
            col(j) = find(rM(j,:),1,"first");
        end
        r = vpa(rM(1:rank(rM),setdiff(1:3,col)));
        if opt.verbose
            fprintf("R%s\n",join(string(col),""))
        end
        switch join(string(col),"")
            case "1"
                r_root = [];
                for eq = subs(equations,[p_var;lin_vars(1:2)],[p_root;-r(1).*lin_vars(2)-r(2);lin_vars(2)])'
                    r_root = unique([r_root,solve(eq,lin_vars(2))]);
                end
                if opt.show_lines
                    p = plot3(repmat(p_root,size(-5:0.1:5)),-r(1).*(-5:0.1:5)-r(2),-5:0.1:5);
                    set(p,"DisplayName","solution","Color","black","LineWidth",1.5)
                    hold on
                end
                q_root = -r(1).*r_root-r(2);
            case "2"
                q_root = [];
                for eq = subs(equations,[p_var;lin_vars(1:2)],[p_root;lin_vars(1);-r(2)])'
                    q_root = unique([q_root,solve(eq,lin_vars(1))]);
                end
                r_root = repmat(-r(2),size(q_root));
                if opt.show_lines
                    p = plot3(repmat(p_root,size(-5:0.1:5)),-5:0.1:5,repmat(-r(2),size(-5:0.1:5)));
                    set(p,"DisplayName","solution","Color","black","LineWidth",1.5)
                    hold on
                end
        end
        % Add all solutions to result vector
        if ~isempty(q_root)
            result = [repmat(p_root,size(q_root)),q_root,r_root];
        end
        % If rank(M)=2 solve system with SVD for numerical stability
    case 2
        [~,S,V] = svd(M);
        [~,index] = min(diag(S));
        result = [p_root,V(1,index)/V(3,index),V(2,index)/V(3,index)];
end
end