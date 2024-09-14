function [v_sol,w_sol] = rank5_4C3_12345(r)
%RANK5_12345 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (5,1) {mustBeReal}
end

if abs(r(1,1) + r(4,1)^2) < 10^-10 && ...
        abs(r(3,1) + r(5,1)^2) < 10^-10 &&...
        abs(r(2,1) + r(4,1)*r(5,1)) < 10^-10
    v_sol = -r(4,1);
    w_sol = -r(5,1);
else
    v_sol = [];
    w_sol = [];
end
end

