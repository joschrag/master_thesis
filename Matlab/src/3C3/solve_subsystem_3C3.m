function [result] = solve_subsystem_3C3(M,p_root,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
    M (6,6) sym;
    p_root (1,1) sym;
    opt.plot_subspace {mustBeInRange(opt.plot_subspace,0,1)} = 1;
    opt.tolerance {mustBeReal,mustBePositive} = 10^-10;
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 1;
    opt.complex (1,1) {mustBeNumericOrLogical} = false;
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
    col(j) = find(abs(rM(j,:)-1) < opt.tolerance,1,"first");
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
        [u_root,v_root] = rank3_3C3_12(r,opt.complex);
    case "13"
        [u_root,v_root] = rank3_3C3_13(r,opt.complex);
    case "14"
        [u_root,v_root] = rank3_3C3_14(r,opt.complex);
    case "15"
        [u_root,v_root] = rank3_3C3_15(r,opt.complex);
    case "23"
        [u_root,v_root] = rank3_3C3_23(r,opt.complex);
    case "24"
        [u_root,v_root] = rank3_3C3_24(r,opt.complex);
    case "25"
        [u_root,v_root] = rank3_3C3_25(r,opt.complex);
    case "34"
        [u_root,v_root] = rank3_3C3_34(r,opt.complex);
    case "35"
        [u_root,v_root] = rank3_3C3_35(r);
    case "45"
        [u_root,v_root] = rank3_3C3_45(r);
    case "123"
        [u_root,v_root] = rank3_3C3_123(r,opt.complex);
    case "124"
        [u_root,v_root] = rank3_3C3_124(r,opt.complex);
    case "125"
        [u_root,v_root] = rank3_3C3_125(r,opt.complex);
    case "134"
        [u_root,v_root] = rank3_3C3_134(r,opt.complex);
    case "135"
        [u_root,v_root] = rank3_3C3_135(r,opt.complex);
    case "145"
        [u_root,v_root] = rank3_3C3_145(r);
    case "234"
        [u_root,v_root] = rank3_3C3_234(r,opt.complex);
    case "235"
        [u_root,v_root] = rank3_3C3_235(r,opt.complex);
    case "245"
        [u_root,v_root] = rank3_3C3_245(r);
    case "345"
        [u_root,v_root] = rank3_3C3_345(r);
    case "1234"
        [u_root,v_root] = rank4_3C3_1234(r,opt.complex);
    case "1245"
        [u_root,v_root] = rank4_3C3_1245(r);
    case "1345"
        [u_root,v_root] = rank4_3C3_1345(r);
    case "2345"
        [u_root,v_root] = rank4_3C3_2345(r);
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