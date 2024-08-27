function params = hyperboloid1_nf(eig_vals,P,offsets,const,opt)
%ELLIPSOID_NF Summary of this function goes here
%   Detailed explanation goes here
arguments
    eig_vals (1,3) {mustBeReal};
    P (3,3) {mustBeReal};
    offsets (3,1) {mustBeReal};
    const (1,1) {mustBeReal};
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
end
if const > 0
    const = -const;
    eig_vals = -eig_vals;
end
n_idx_p = (1:3).*(eig_vals>0);
n_idx_m = (1:3).*(eig_vals<0);
n_idx = [n_idx_p(n_idx_p~=0),n_idx_m(n_idx_m~=0)];
eig_vals = eig_vals(n_idx);
tmp = 1./sqrt(abs(eig_vals)./-const);
P = P(:,n_idx);
offsets = offsets(n_idx);
params = get_hyperboloid_1_params(tmp(1),tmp(2),tmp(3),P,offsets,plotRange=opt.plotRange);
end