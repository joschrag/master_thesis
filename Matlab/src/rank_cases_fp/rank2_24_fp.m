function [u_sol,v_sol] = rank2_24_fp(r,prime)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (2,3) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePositive,mustBeInteger}
end
if r(1,2) == 0
    if FF(r(1,3),prime) == -FF(r(2,3),prime)^2
        u_0 = FF(-r(2,3),prime).value;
        t_0 = 1:prime;
        tmp = zeros(numel(t_0),2);
        for i=1:numel(t_0)
            tmp(i,:) = [t_0(i),u_0];
        end
        u_sol = [tmp(:,1)];
        v_sol = [tmp(:,2)];
    else
        v_sol = [];
        u_sol = [];
    end
else
    v_sol = FF(-r(2,3),prime).value;
    t_0 = -FF(r(1,3)+r(2,3)^2,prime)*FF(r(1,2),prime)^(-1);
    u_sol = t_0.value;
end
end