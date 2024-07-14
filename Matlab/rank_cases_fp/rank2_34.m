function [u_sol,v_sol] = rank2_34(r,prime)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (2,3) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePositive,mustBeInteger}
end
u_sol = FF(-r(1,3),prime).value;
v_sol = FF(-r(2,3),prime).value;
end