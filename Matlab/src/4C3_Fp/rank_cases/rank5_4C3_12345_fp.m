function [v_sol,w_sol] = rank5_4C3_12345_fp(r,prime)
%RANK5_12345 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (5,1) {mustBeReal}
    prime (1,1) {mustBeInteger,mustBePositive}
end

if FF(r(1,1) + r(4,1)^2,prime).value == 0 && ...
        FF(r(2,1) + r(4,1)*r(5,1),prime).value == 0 && ...
        FF(r(3,1) + r(5,1)^2,prime).value == 0
    v_sol = FF(-r(4,1),prime).value;
    w_sol =  FF(-r(5,1),prime).value;
else
    v_sol = [];
    w_sol = [];
end
end

