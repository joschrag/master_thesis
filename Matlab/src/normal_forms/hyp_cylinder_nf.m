function params = hyp_cylinder_nf(eig_vals,P,offsets,const)
%HYP_CYLINDER_NF Summary of this function goes here
%   Detailed explanation goes here
n_idx_p = (1:3).*(eig_vals>0)';
n_idx_m = (1:3).*(eig_vals<0)';
n_idx_z = (1:3).*(eig_vals==0)';
n_idx = [n_idx_p(n_idx_p~=0),n_idx_m(n_idx_m~=0),n_idx_z(n_idx_z~=0)];
eig_vals = eig_vals(n_idx);
P = P(:,n_idx);
offsets = offsets(n_idx);
if const > 0
    const = -const;
    eig_vals(1:2) = -eig_vals(2:-1:1);
    P = -P(:,[2,1,3]);
end

tmp = 1./sqrt(abs(eig_vals)./-const);
params = get_hyp_cylinder_params(tmp(1),tmp(2),P,offsets);
end

