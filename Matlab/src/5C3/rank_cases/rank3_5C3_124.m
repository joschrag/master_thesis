function [u_sol,v_sol] = rank3_5C3_124(r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (3,2) {mustBeReal}
end

v_0 = -r(3,2);
u_0 = roots([-1,-r(1,1),-r(1,2)]);
u_0 = u_0(imag(u_0)==0);
control = abs(-v_0^2-r(2,1).*u_0-r(2,2)) < 10^-10;
if ~any(control)
    v_sol = [];
    u_sol = [];
    return
end
u_0 = u_0(control);
tmp = zeros(numel(u_0),2);
for i=1:numel(u_0)
    tmp(i,:) = [u_0(i),v_0];
end
u_sol = tmp(:,1);
v_sol = tmp(:,2);
end

