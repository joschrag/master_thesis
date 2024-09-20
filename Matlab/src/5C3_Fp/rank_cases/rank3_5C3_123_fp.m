function [u_sol,v_sol] = rank3_5C3_123_fp(r,prime)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (3,2) {mustBeReal}
    prime (1,1) {mustBePositive,mustBeInteger}
end

v_sol = get_gf_root([-1,-r(2,:)],prime);
if isempty(v_sol)
    u_sol = [];
    return
end
u_sol = FF(-(r(3,1).*v_sol+r(3,2)),prime).value;

control = FF(-(r(1,1).*v_sol+u_sol.^2+r(1,2)),prime).value == 0;

v_sol = v_sol(control);
u_sol = u_sol(control);
end

