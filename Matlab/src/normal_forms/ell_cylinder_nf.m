function params = ell_cylinder_nf(eig_vals,V,offsets,const)
%ELL_CYLINDER_NF Summary of this function goes here
%   Detailed explanation goes here
n_idx_p = (1:3).*(eig_vals>0)';
n_idx_m = (1:3).*(eig_vals<0)';
n_idx_z = (1:3).*(eig_vals==0)';
n_idx = [n_idx_p(n_idx_p~=0),n_idx_m(n_idx_m~=0),n_idx_z(n_idx_z~=0)];
eig_vals = eig_vals(n_idx);
V = V(:,n_idx);

offsets = offsets(n_idx);
if const > 0
    const = -const;
    eig_vals = -eig_vals;
end
if sign(eig_vals(1)) == sign(const)
    error("Imaginary elliptic Cylinder")
end
tmp = 1./sqrt(eig_vals./-const);
params = get_ell_cylinder_params(tmp(1),tmp(2),V,offsets);
end

