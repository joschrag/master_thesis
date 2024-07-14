function params = ell_cone_nf(eig_vals,V,offsets)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n_idx_p = (1:3).*(eig_vals>0)';
n_idx_m = (1:3).*(eig_vals<0)';
n_idx = [n_idx_p(n_idx_p~=0),n_idx_m(n_idx_m~=0)];
factor = 1;
if sum(n_idx_m~=0) == 2
    n_idx = [n_idx_m(n_idx_m~=0),n_idx_p(n_idx_p~=0)];
    factor = -1;
end
eig_vals = factor.*eig_vals(n_idx);
tmp = 1./sqrt(abs(eig_vals));
V = factor.*V(:,n_idx);
offsets = factor.*offsets(n_idx);
params = get_ell_cone_params(tmp(1),tmp(2),tmp(3),V,offsets);
end

