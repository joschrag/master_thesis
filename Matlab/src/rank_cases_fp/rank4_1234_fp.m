function [u_sol,v_sol] = rank4_1234_fp(r,prime)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (4,1) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePositive,mustBeInteger}
end


if (FF(-r(1,1),prime) == FF(r(3,1),prime)^2) && ...
        (FF(-r(2,1),prime) == FF(r(4,1),prime)^2)
    u_sol = mod(-r(3,1),prime);
    v_sol = mod(-r(4,1),prime);
else
    u_sol = [];
    v_sol = [];
end
end