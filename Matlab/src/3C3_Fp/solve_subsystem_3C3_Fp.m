function [result] = solve_subsystem_3C3_Fp(M,p_root,prime,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
    M (1,1) FF;
    p_root (1,1) sym {mustBeReal};
    prime (1,1) {mustBeInteger,mustBePositive};
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 1;
end
result = [];
if rank(M)==0
    warning("Zero matrix for p_root %i",p_root)
    return
end
if rank(M) == 6
    warning("Precision of root %i is too low!",p_root)
    return
end
rM = rref(M);
if rank(M) ~= rank(rM)
    warning("Numerical error calculating rref matrix for p_root %i",p_root)
    return
end
col = zeros(1,rank(rM));
for j=1:rank(rM)
    col(j) = find(rM(j,:),1,"first");
end
m = size(M.value,1);
vec = 1:m;
if rank(rM) ~= rank(M)
    warning("Numerical error!")
    return
end
r = rM(1:rank(rM),setdiff(vec,col));
if opt.verbose
    fprintf("R%s\n",join(string(col),""))
end
switch join(string(col),"")
    case "1"
        [u_root,v_root] = rank1_3C3_1_Fp(r,prime);
    case "2"
        [u_root,v_root] = rank1_3C3_2_Fp(r,prime);
    case "3"
        [u_root,v_root] = rank1_3C3_3_Fp(r,prime);
    case "4"
        [u_root,v_root] = rank1_3C3_4_Fp(r,prime);
    case "12"
        [u_root,v_root] = rank2_3C3_12_Fp(r,prime);
    case "13"
        [u_root,v_root] = rank2_3C3_13_Fp(r,prime);
    case "14"
        [u_root,v_root] = rank2_3C3_14_Fp(r,prime);
    case "23"
        [u_root,v_root] = rank2_3C3_23_Fp(r,prime);
    case "24"
        [u_root,v_root] = rank2_3C3_24_Fp(r,prime);
    case "34"
        [u_root,v_root] = rank2_3C3_34_Fp(r,prime);
    case "123"
        [u_root,v_root] = rank3_3C3_123_Fp(r,prime);
    case "124"
        [u_root,v_root] = rank3_3C3_124_Fp(r,prime);
    case "134"
        [u_root,v_root] = rank3_3C3_134_Fp(r,prime);
    case "234"
        [u_root,v_root] = rank3_3C3_234_Fp(r,prime);
    case "1234"
        [u_root,v_root] = rank4_3C3_1234_Fp(r,prime);
    case "1245"
        [u_root,v_root] = rank4_3C3_1234_Fp(r,prime);
    case "1345"
        [u_root,v_root] = rank4_3C3_1234_Fp(r,prime);
    case "2345"
        [u_root,v_root] = rank4_3C3_1234_Fp(r,prime);
    case "12345"
        [u_root,v_root] = rank5_3C3_12345_Fp(r,prime);
    otherwise
        return
        % fprintf("System of equations is singular.\n");
end
% Add solutions to result vector
if ~isempty(u_root)
    result=[repmat(p_root,numel(u_root),1),reshape(u_root,[],1),reshape(v_root,[],1)];
end
end