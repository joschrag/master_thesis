function [u_sol,v_sol] = rank3_5C3_134(r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (3,2) {mustBeReal}
end

if abs(-r(1,1)*r(3,2)^2-r(1,2)-r(2,2)^2) < 10^-10
    v_sol = -r(3,2);
    u_sol = -r(2,2);
else
    v_sol = [];
    u_sol = [];
    return
end
end

