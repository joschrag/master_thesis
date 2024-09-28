function result = solve_subsystem_E3Q3_Fp(M,p_root,prime)
%SOLVE_SUBSYSTEM_E3Q3_FP Solve the resulting system of equations for all cases.
%   For all cases of rref(M) this function solves the system of equations and
%   returns the solution.
arguments
    M FF {mustBeSizeFF(M,[3,3])};
    p_root (1,1) {mustBeInteger};
    prime (1,1) {mustBePrime};
end
result = [];
if isequaln(M.value,zeros(3))
    return
end

rM = rref(M);
% Construct array of indepedent columns
m = zeros(1,rank(rM));
for j=1:rank(rM)
    m(j) = find(rM(j,:),1,"first");
end
r = vpa(rM(1:rank(rM),setdiff(1:3,m)));
% Solve each of the cases accordingly
switch join(string(m),"")
    case "1"
        for eq = subs(equations,[p_var;lin_vars(1:2)],[p_root;-r(1).*lin_vars(2)-r(2);lin_vars(2)])'
            r_root = get_gf_root(coeffs(eq,lin_vars.value(2),"All"),prime);
        end
        q_root = FF(-r(1).*r_root-r(2),prime).value;
    case "2"
        for j=1:3
            eq = subs(equations(j),[p_var;lin_vars.value(1:2)],[p_root;lin_vars.value(1);-r(2)]);
            q_root = get_gf_root(coeffs(eq,lin_vars.value(1),"All"),prime);
        end
        r_root = repmat(FF(-r(2),prime).value,size(q_root));
    case "12"
        o_sols = FF(-r',prime).value;
        q_root = o_sols(1);
        r_root = o_sols(2);
end
% Add solutions to result vector
if ~isempty(q_root)
    result = [repmat(p_root,size(q_root)),q_root,r_root];
end
end