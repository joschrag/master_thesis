function params = line_nf(eig_vals,P,offsets)
%LINE_NF Summary of this function goes here
%   Detailed explanation goes here
n_idx_p = (1:3).*(eig_vals>0)';
n_idx_m = (1:3).*(eig_vals<0)';
n_idx_z = (1:3).*(eig_vals==0)';
n_idx = [n_idx_p(n_idx_p~=0),n_idx_m(n_idx_m~=0),n_idx_z(n_idx_z~=0)];
P = P(:,n_idx);
offsets = offsets(n_idx);
params = get_line_params(P,offsets);
end

