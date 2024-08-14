function params = par_cylinder_nf(eig_vals,P,offsets,d0,d1,const)
%PAR_CYLINDER_NF Summary of this function goes here
%   Detailed explanation goes here
n_idx_p = (1:3).*(eig_vals>0)';
n_idx_m = (1:3).*(eig_vals<0)';
n_idx_z = (1:3).*(eig_vals==0)';
n_idx = [n_idx_p(n_idx_p~=0),n_idx_m(n_idx_m~=0),n_idx_z(n_idx_z~=0)];
eig_vals = eig_vals(n_idx);
P = P(:,n_idx);
offsets = offsets(n_idx);
if d1 > 0
    d0 = -d0;
    d1 = -d1;
    const = -const;
    eig_vals = -eig_vals;
    d = d1;
elseif d1 < 0
        d = d1;
elseif d1 == 0
    if d0 > 0
    d0 = -d0;
    const = -const;
    eig_vals = -eig_vals;
    d = d0;
    else
    d = d0;
    end
end
offsets(3) = -const/d;
tmp = 1./sqrt(eig_vals(1)./d);
params = get_par_cylinder_params(tmp(1),d0/d,P,offsets);
end

