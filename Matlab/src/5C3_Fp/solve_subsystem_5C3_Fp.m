function result = solve_subsystem_5C3_Fp(M,p_root,prime,opt)
%SOLVE_SUBSYSTEM_5C3_FP Solve the resulting system of equations for all cases.
%   For all cases of rref(M) this function solves the system of equations and
%   returns the solution.
arguments
    M FF {mustBeSizeFF(M,[5,5])};
    p_root (1,1) {mustBeInteger};
    prime (1,1) {mustBePrime};
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 1;
end
result = [];
if rank(M)==0
    warning("Zero matrix for p_root %i",p_root)
    return
end
if rank(M) == 5
    warning("Precision of root %i is too low!",p_root)
    return
end
rM = rref(M);
if rank(M) ~= rank(rM)
    warning("Numerical error calculating rref matrix for p_root %i",p_root)
    return
end
% Construct array of indepedent columns
col = zeros(1,rank(rM));
for j=1:rank(rM)
    col(j) = find(rM(j,:),1,"first");
end
r = vpa(rM(1:rank(rM),setdiff(1:5,col)));

if opt.verbose
    fprintf("R%s\n",join(string(col),""))
end

switch join(string(col),"")
    case "1"
        [q_root,r_root] = rank1_5C3_1_Fp(r,prime);
    case "2"
        [q_root,r_root] = rank1_5C3_2_Fp(r,prime);
    case "3"
        [q_root,r_root] = rank1_5C3_3_Fp(r,prime);
    case "4"
        [q_root,r_root] = rank1_5C3_4_Fp(r,prime);
    case "12"
        [q_root,r_root] = rank2_5C3_12_Fp(r,prime);
    case "13"
        [q_root,r_root] = rank2_5C3_13_Fp(r,prime);
    case "14"
        [q_root,r_root] = rank2_5C3_14_Fp(r,prime);
    case "23"
        [q_root,r_root] = rank2_5C3_23_Fp(r,prime);
    case "24"
        [q_root,r_root] = rank2_5C3_24_Fp(r,prime);
    case "34"
        [q_root,r_root] = rank2_5C3_34_Fp(r,prime);
    case "123"
        [q_root,r_root] = rank3_5C3_123_Fp(r,prime);
    case "124"
        [q_root,r_root] = rank3_5C3_124_Fp(r,prime);
    case "134"
        [q_root,r_root] = rank3_5C3_134_Fp(r,prime);
    case "234"
        [q_root,r_root] = rank3_5C3_234_Fp(r,prime);
    case "1234"
        [q_root,r_root] = rank4_5C3_1234_Fp(r,prime);
    otherwise
        fprintf("System of equations is singular.\n");
end
% Add all solutions to result vector
if ~isempty(q_root)
    result = [repmat(p_root,size(q_root)),q_root,r_root];
end
end