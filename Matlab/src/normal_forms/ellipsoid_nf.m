function params = ellipsoid_nf(eig_vals,P,offsets,const)
%ELLIPSOID_NF Summary of this function goes here
%   Detailed explanation goes here
arguments
    eig_vals (1,3) {mustBeReal};
    P (3,3) {mustBeReal};
    offsets (3,1) {mustBeReal};
    const (1,1) {mustBeReal};
end
if const > 0
    const = -const;
    eig_vals = -eig_vals;
end
tmp = 1./sqrt(eig_vals./-const);
params = get_ellipsoid_params(tmp(1),tmp(2),tmp(3),P,offsets);
end

