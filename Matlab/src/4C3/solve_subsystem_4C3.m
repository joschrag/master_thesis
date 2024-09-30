function [result] = solve_subsystem_4C3(M,p_root,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
    M (6,6) sym {mustBeReal};
    p_root (1,1) sym {mustBeReal};
    opt.plot_subspace {mustBeInRange(opt.plot_subspace,0,1)} = 1;
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 1;
    opt.tolerance (1,1) {mustBeReal,mustBePositive} = 10^-10;
end
result = [];
if rank(M)==0
    warning("Zero matrix for p_root %i",p_root)
    return
end
if rank(double(M),opt.tolerance) == 6
    warning("Precision of root %f is too low!",p_root)
    return
end
rM = rref(double(M), 1e-6);
if rank(M) ~= rank(rM)
    warning("Numerical error calculating rref matrix for p_root %i",p_root)
    return
end
col = zeros(1,rank(rM));
for j=1:rank(rM)
    col(j) = find(rM(j,:),1,"first");
end
m = size(M,1);
vec = 1:m;
r = vpa(rM(1:rank(rM),setdiff(vec,col)));
if opt.verbose
    fprintf("R%s\n",join(string(col),""))
end
switch join(string(col),"")
    case "1"
        [u_root,v_root] = rank1_4C3_1(r);
    case "2"
        [u_root,v_root] = rank1_4C3_2(r);
    case "3"
        [u_root,v_root] = rank1_4C3_3(r);
    case "4"
        [u_root,v_root] = rank1_4C3_4(r);
    case "12"
        [u_root,v_root] = rank2_4C3_12(r);
    case "13"
        [u_root,v_root] = rank2_4C3_13(r);
    case "14"
        [u_root,v_root] = rank2_4C3_14(r);
    case "23"
        [u_root,v_root] = rank2_4C3_23(r);
    case "24"
        [u_root,v_root] = rank2_4C3_24(r);
    case "34"
        [u_root,v_root] = rank2_4C3_34(r);
    case "123"
        [u_root,v_root] = rank3_4C3_123(r);
    case "124"
        [u_root,v_root] = rank3_4C3_124(r);
    case "134"
        [u_root,v_root] = rank3_4C3_134(r);
    case "234"
        [u_root,v_root] = rank3_4C3_234(r);
    case "1234"
        [u_root,v_root] = rank4_4C3_1234(r);
    case "1245"
        [u_root,v_root] = rank4_4C3_1234(r);
    case "1345"
        [u_root,v_root] = rank4_4C3_1234(r);
    case "2345"
        [u_root,v_root] = rank4_4C3_1234(r);
    case "12345"
        [u_root,v_root] = rank5_4C3_12345(r);
    otherwise
        return
        % fprintf("System of equations is singular.\n");
end
% Add solutions to result vector
if ~isempty(u_root)
    result=[repmat(p_root,numel(u_root),1),reshape(u_root,[],1),reshape(v_root,[],1)];
end
end