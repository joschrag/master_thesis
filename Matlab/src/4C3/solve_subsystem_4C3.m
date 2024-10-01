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
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "2"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "3"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "4"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "12"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "13"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "14"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "23"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "24"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "34"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "123"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "124"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "134"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "234"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "235"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "245"
        warning("Function to solve R%s not implemented!",join(string(col),""))
    case "345"
        warning("Function to solve R%s not implemented!",join(string(col),""))
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