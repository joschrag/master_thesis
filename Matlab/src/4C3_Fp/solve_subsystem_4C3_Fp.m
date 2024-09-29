function result = solve_subsystem_4C3_Fp(M,p_root,prime)
%SOLVE_SUBSYSTEM_4C3_FP Solve the resulting system of equations for all cases.
%   For all cases of rref(M) this function solves the system of equations and
%   returns the solution.
arguments
    M FF {mustBeSizeFF(M,[6,6])};
    p_root (1,1) {mustBeInteger};
    prime (1,1) {mustBePrime};
end
result = [];
if isequaln(M.value,zeros(6))
    return
end
rM = rref(M);
% Construct array of indepedent columns
col = zeros(1,rank(rM));
for j=1:rank(rM)
    col(j) = find(rM(j,:),1,"first");
end
r = vpa(rM(1:rank(rM),setdiff(1:6,col)));

switch join(string(col),"")
    case "1"
        [u_root,v_root] = rank1_4C3_1_fp(r,prime);
    case "2"
        [u_root,v_root] = rank1_4C3_2_fp(r,prime);
    case "3"
        [u_root,v_root] = rank1_4C3_3_fp(r,prime);
    case "4"
        [u_root,v_root] = rank1_4C3_4_fp(r,prime);
    case "12"
        [u_root,v_root] = rank2_4C3_12_fp(r,prime);
    case "13"
        [u_root,v_root] = rank2_4C3_13_fp(r,prime);
    case "14"
        [u_root,v_root] = rank2_4C3_14_fp(r,prime);
    case "23"
        [u_root,v_root] = rank2_4C3_23_fp(r,prime);
    case "24"
        [u_root,v_root] = rank2_4C3_24_fp(r,prime);
    case "34"
        [u_root,v_root] = rank2_4C3_34_fp(r,prime);
    case "123"
        [u_root,v_root] = rank3_4C3_123_fp(r,prime);
    case "124"
        [u_root,v_root] = rank3_4C3_124_fp(r,prime);
    case "134"
        [u_root,v_root] = rank3_4C3_134_fp(r,prime);
    case "234"
        [u_root,v_root] = rank3_4C3_234_fp(r,prime);
    case "1234"
        [u_root,v_root] = rank4_4C3_1234_fp(r,prime);
    case "1245"
        [u_root,v_root] = rank4_4C3_1234_fp(r,prime);
    case "1345"
        [u_root,v_root] = rank4_4C3_1234_fp(r,prime);
    case "2345"
        [u_root,v_root] = rank4_4C3_1234_fp(r,prime);
    case "12345"
        [u_root,v_root] = rank5_4C3_12345_fp(r,prime);
    otherwise
        fprintf("System of equations is singular.\n");
        return
end
if ~isempty(u_root)
    result = [repmat(p_root,size(u_root)),u_root,v_root];
end
end