function params = plane_nf(eig_vals,P,offsets,opt)
%PLANE_NF Summary of this function goes here
%   Detailed explanation goes here
arguments
    eig_vals (1,3) {mustBeReal};
    P (3,3) {mustBeReal};
    offsets (3,1) {mustBeReal};
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
    opt.tolerance (1,1) {mustBeReal,mustBePositive} = 10^-10;
end
n_idx_p = (1:3).*(eig_vals>opt.tolerance)';
n_idx_m = (1:3).*(eig_vals<-opt.tolerance)';
n_idx_z = (1:3).*(abs(eig_vals)<opt.tolerance)';
n_idx = [n_idx_p(n_idx_p~=0),n_idx_m(n_idx_m~=0),n_idx_z(n_idx_z~=0)];
P = P(:,n_idx);
offsets = offsets(n_idx);
params = get_plane_params(P,offsets,plotRange=opt.plotRange);
end

