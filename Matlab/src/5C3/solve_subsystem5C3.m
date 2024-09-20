function [result] = solve_subsystem5C3(M,p_root,idx,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
    M (5,5) sym {mustBeReal};
    p_root (1,1) sym {mustBeReal};
    idx (1,3) {mustBeInteger};
    opt.plot_subspace {mustBeInRange(opt.plot_subspace,0,1)} = 1;
end
result = [];
rM = rref(double(M), 1e-6);
col = zeros(1,rank(rM));
for j=1:rank(rM)
    col(j) = find(rM(j,:),1,"first");
end
r = vpa(rM(1:rank(rM),setdiff(1:5,col)));
output = zeros(1, 5);
vec = 1:5;
output(col) = vec(1:rank(rM));

remaining_elements = vec(rank(rM) + 1:end);
remaining_indices = setdiff(1:5, col);
output(remaining_indices) = remaining_elements;
base = [-r;eye(5-rank(rM))];
base = base(output,:);
switch join(string(col),"")
    case "1"
        [u_root,v_root] = rank1_5C3_1(r);
    case "2"
        [u_root,v_root] = rank1_5C3_2(r);
    case "3"
        [u_root,v_root] = rank1_5C3_3(r);
    case "4"
        [u_root,v_root] = rank1_5C3_4(r);
    case "12"
        [u_root,v_root] = rank2_5C3_12(r);
    case "13"
        [u_root,v_root] = rank2_5C3_13(r);
    case "14"
        [u_root,v_root] = rank2_5C3_14(r);
    case "23"
        [u_root,v_root] = rank2_5C3_23(r);
    case "24"
        [u_root,v_root] = rank2_5C3_24(r);
    case "34"
        [u_root,v_root] = rank2_5C3_34(r);
    case "123"
        [u_root,v_root] = rank3_5C3_123(r);
    case "124"
        [u_root,v_root] = rank3_5C3_124(r);
    case "134"
        [u_root,v_root] = rank3_5C3_134(r);
    case "234"
        [u_root,v_root] = rank3_5C3_234(r);
    case "1234"
        [u_root,v_root] = rank4_5C3_1234(r);
    otherwise
        return
        % fprintf("System of equations is singular.\n");
end
if ~isempty(u_root)
    result=[repmat(p_root,numel(u_root),1),reshape(u_root,[],1),reshape(v_root,[],1)];
    result = result(:,idx);
end
if opt.plot_subspace
    plot_subspace(u_root,v_root,rank(rM),base)
end
end