function params = hyperboloid1_nf(eig_vals,V,offsets,const)
%ELLIPSOID_NF Summary of this function goes here
%   Detailed explanation goes here
if const > 0
    const = -const;
    eig_vals = -eig_vals;
end
n_idx_p = (1:3).*(eig_vals>0)';
n_idx_m = (1:3).*(eig_vals<0)';
n_idx = [n_idx_p(n_idx_p~=0),n_idx_m(n_idx_m~=0)];
eig_vals = eig_vals(n_idx);
tmp = 1./sqrt(abs(eig_vals)./-const);
V = V(:,n_idx);
offsets = offsets(n_idx);
params = get_hyperboloid_1_params(tmp(1),tmp(2),tmp(3),V,offsets);
end