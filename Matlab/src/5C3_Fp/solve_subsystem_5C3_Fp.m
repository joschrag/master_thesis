function result = solve_subsystem_5C3_Fp(M,p_root,prime)
%SOLVE_SUBSYSTEM_5C3_FP Solve the resulting system of equations for all cases.
%   For all cases of rref(M) this function solves the system of equations and
%   returns the solution.
arguments
    M FF {mustBeSizeFF(M,[5,5])};
    p_root (1,1) {mustBeInteger};
    prime (1,1) {mustBePrime};
end
if isequaln(M.value,zeros(5))
    return
end
result = [];
rM = rref(M);
% Construct array of indepedent columns
m = zeros(1,rank(rM));
for j=1:rank(rM)
    m(j) = find(rM(j,:),1,"first");
end
r = vpa(rM(1:rank(rM),setdiff(1:5,m)));
switch join(string(m),"")
    case "1"
        [q_root,r_root] = rank1_5C3_1_fp(r,prime);
    case "2"
        [q_root,r_root] = rank1_5C3_2_fp(r,prime);
    case "3"
        [q_root,r_root] = rank1_5C3_3_fp(r,prime);
    case "4"
        [q_root,r_root] = rank1_5C3_4_fp(r,prime);
    case "12"
        [q_root,r_root] = rank2_5C3_12_fp(r,prime);
    case "13"
        [q_root,r_root] = rank2_5C3_13_fp(r,prime);
    case "14"
        [q_root,r_root] = rank2_5C3_14_fp(r,prime);
    case "23"
        [q_root,r_root] = rank2_5C3_23_fp(r,prime);
    case "24"
        [q_root,r_root] = rank2_5C3_24_fp(r,prime);
    case "34"
        [q_root,r_root] = rank2_5C3_34_fp(r,prime);
    case "123"
        [q_root,r_root] = rank3_5C3_123_fp(r,prime);
    case "124"
        [q_root,r_root] = rank3_5C3_124_fp(r,prime);
    case "134"
        [q_root,r_root] = rank3_5C3_134_fp(r,prime);
    case "234"
        [q_root,r_root] = rank3_5C3_234_fp(r,prime);
    case "1234"
        [q_root,r_root] = rank4_5C3_1234_fp(r,prime);
    otherwise
        fprintf("System of equations is singular.\n");
end
% Add all solutions to result vector
if ~isempty(q_root)
    result = [repmat(p_root,size(q_root)),q_root,r_root];
end
end