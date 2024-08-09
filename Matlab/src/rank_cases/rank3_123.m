function [u_sol,v_sol] = rank3_123(r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (3,2) {mustBeReal}
end

v_sol = roots([-1,-r(2,:)]);
v_sol = unique(v_sol(imag(v_sol)==0));
if isempty(v_sol)
    u_sol = [];
    return
end
u_sol = -(r(3,1).*v_sol+r(3,2));
control = abs(-(r(1,1).*v_sol+u_sol.^2+r(1,2))) < 10^-10;

v_sol = v_sol(control);
u_sol = u_sol(control);
end

