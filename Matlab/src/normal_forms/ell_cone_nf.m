function params = ell_cone_nf(eig_vals,P,offsets,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
    eig_vals (1,3) {mustBeReal};
    P (3,3) {mustBeReal};
    offsets (3,1) {mustBeReal};
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
end
n_idx_p = (1:3).*(eig_vals>0);
n_idx_m = (1:3).*(eig_vals<0);
n_idx = [n_idx_p(n_idx_p~=0),n_idx_m(n_idx_m~=0)];
factor = 1;
if sum(n_idx_m~=0) == 2
    n_idx = [n_idx_m(n_idx_m~=0),n_idx_p(n_idx_p~=0)];
    factor = -1;
end
eig_vals = factor.*eig_vals(n_idx);
tmp = 1./sqrt(abs(eig_vals));
P = factor.*P(:,n_idx);
offsets = factor.*offsets(n_idx);
params = get_ell_cone_params(tmp(1),tmp(2),tmp(3),P,offsets,plotRange=opt.plotRange);
end

