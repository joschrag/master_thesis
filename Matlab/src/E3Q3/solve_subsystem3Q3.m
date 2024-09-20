function [result] = solve_subsystem3Q3(M,p_root,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
    M (3,3) sym {mustBeReal};
    p_root (1,1) sym {mustBeReal};
    opt.show_lines (1,1) {mustBeNumericOrLogical} = false;
    opt.tolerance (1,1) {mustBeReal,mustBePositive} = 10^-10;
end
switch rank(double(M),opt.tolerance)
    case 1
        rM = rref(double(M), opt.tolerance);
        col = zeros(1,rank(rM));
        for j=1:rank(rM)
            col(j) = find(rM(j,:),1,"first");
        end
        r = vpa(rM(1:rank(rM),setdiff(1:3,col)));
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
        if ~isempty(q_root)
            for root=[q_root,r_root]'
                result=[p_root,root'];
            end
        end
    case 2
        [~,S,V] = svd(M);
        [~,index] = min(diag(S));
        result = [p_root,V(1,index)/V(3,index),V(2,index)/V(3,index)];
end
end