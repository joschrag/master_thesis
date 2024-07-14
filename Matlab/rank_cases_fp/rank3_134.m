function [u_sol,v_sol] = rank3_134(r,prime)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (3,2) {mustBeReal}
    prime (1,1) {mustBePositive,mustBeInteger}
end

if FF(-r(1,1)*r(3,2)^2-r(1,2)-r(2,2)^2,prime).value == 0
    v_sol = FF(-r(3,2),prime).value;
    u_sol = FF(-r(2,2),prime).value;
else
    v_sol = [];
    u_sol = [];
    return
end
end

