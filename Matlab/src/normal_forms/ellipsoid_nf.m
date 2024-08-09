function params = ellipsoid_nf(eig_vals,V,offsets,const)
%ELLIPSOID_NF Summary of this function goes here
%   Detailed explanation goes here
if const > 0
    const = -const;
    eig_vals = -eig_vals;
end
tmp = 1./sqrt(eig_vals./-const);
params = get_ellipsoid_params(tmp(1),tmp(2),tmp(3),V,offsets);
end

