function [u_sol,v_sol] = rank4_1234(r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (4,1) {mustBeReal}
end


if abs(r(1,1) + r(3,1)^2) < 10^-10 && ...
        abs(r(2,1) + r(4,1)^2) < 10^-10
    u_sol = -r(3,1);
    v_sol = -r(4,1);
else
    u_sol = [];
    v_sol = [];
end
end