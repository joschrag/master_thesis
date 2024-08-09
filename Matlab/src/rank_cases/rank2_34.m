function [u_sol,v_sol] = rank2_34(r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (2,3) {mustBeReal}
end
u_sol = -r(1,3);
v_sol = -r(2,3);
end