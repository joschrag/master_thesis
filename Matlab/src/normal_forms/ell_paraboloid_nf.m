function params = ell_paraboloid_nf(eig_vals, V, offsets, d0,const)
n_idx_p = (1:3).*(eig_vals>0)';
n_idx_m = (1:3).*(eig_vals<0)';
n_idx_z = (1:3).*(eig_vals==0)';
n_idx = [n_idx_p(n_idx_p~=0),n_idx_m(n_idx_m~=0),n_idx_z(n_idx_z~=0)];
V = V(:,n_idx);
eig_vals = eig_vals(n_idx);
offsets = offsets(n_idx);
offsets(3) = const/d0;
dir=1;
if sign(eig_vals(1)) == sign(d0)
    dir = -1;
end
tmp = 1./sqrt(abs(eig_vals(eig_vals~=0)./d0));
params = get_ell_paraboloid_params(tmp(1),tmp(2),V,offsets,dir);
end