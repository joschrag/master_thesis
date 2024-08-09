function [u_sol,v_sol] = rank2_23(r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (2,3) {mustBeReal}
end
v_0 = roots([-1,-r(1,2),-r(1,3)]);
v_0 = v_0(imag(v_0)==0);
if ~isempty(v_0)
    v_sol = v_0;
    u_sol = -r(2,2).*v_0-r(2,3);
else
    v_sol = [];
    u_sol = [];
end
end