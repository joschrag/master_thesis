function [result] = solve_subsystem_3C3(M,p_root,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
    M (6,6) sym {mustBeReal};
    p_root (1,1) sym {mustBeReal};
    opt.plot_subspace {mustBeInRange(opt.plot_subspace,0,1)} = 1;
    opt.tolerance {mustBeReal,mustBePositive} = 10^-10;
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 1;
end
result = [];
if rank(M)==0
    warning("Zero matrix for p_root %f",p_root)
    return
end
if rank(double(M),opt.tolerance) == 6
    warning("Precision of root %f is too low!",p_root)
    return
end
rM = rref(double(M), opt.tolerance);
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
if rank(rM) ~= rank(M)
    warning("Numerical error!")
    return
end
r = vpa(rM(1:rank(rM),setdiff(vec,col)));
output = zeros(1, m);
output(col) = vec(1:rank(rM));

remaining_elements = vec(rank(rM) + 1:end);
remaining_indices = setdiff(vec, col);
output(remaining_indices) = remaining_elements;
base = [-r;eye(m-rank(rM))];
base = base(output,:);
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
        [u_root,v_root] = rank3_3C3_123(r);
    case "124"
        [u_root,v_root] = rank3_3C3_124(r);
    case "134"
        [u_root,v_root] = rank3_3C3_134(r);
    case "234"
        [u_root,v_root] = rank3_3C3_234(r);
    case "1234"
        [u_root,v_root] = rank4_3C3_1234(r);
    case "1245"
        [u_root,v_root] = rank4_3C3_1234(r);
    case "1345"
        [u_root,v_root] = rank4_3C3_1234(r);
    case "2345"
        [u_root,v_root] = rank4_3C3_1234(r);
    case "12345"
        [u_root,v_root] = rank5_3C3_12345(r);
    otherwise
        return
        % fprintf("System of equations is singular.\n");
end
% Add solutions to result vector
if ~isempty(u_root)
    result=[repmat(p_root,numel(u_root),1),reshape(u_root,[],1),reshape(v_root,[],1)];
end
% Plot subspace geometry
if opt.plot_subspace
    plot_subspace3C3(u_root,v_root,col,base)
end
end