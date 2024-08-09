function params = par_planes_nf(eig_vals,P,offsets,const)
%PAR_PLANES_NF Summary of this function goes here
%   Detailed explanation goes here
n_idx_p = (1:3).*(eig_vals>0)';
n_idx_m = (1:3).*(eig_vals<0)';
n_idx_z = (1:3).*(eig_vals==0)';
n_idx = [n_idx_p(n_idx_p~=0),n_idx_m(n_idx_m~=0),n_idx_z(n_idx_z~=0)];
eig_vals = eig_vals(n_idx);
P = P(:,n_idx);
offsets = offsets(n_idx);
if const > 0
    eig_vals = -eig_vals;
    const = -const;
end
if sign(eig_vals(1)) == sign(const)
    error("Imaginary parallel planes.")
end
tmp = 1./sqrt(eig_vals(1)./-const);
params = get_par_plane_params(tmp(1),P,offsets);
end

