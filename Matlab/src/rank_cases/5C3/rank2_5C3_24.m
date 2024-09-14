function [u_sol,v_sol] = rank2_5C3_24(r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (2,3) {mustBeReal}
end
if r(1,2) == 0
    if r(1,3) == -r(2,3)^2
        v_0 = -r(2,3);
        u_0 = -10:0.01:10;
        tmp = zeros(numel(u_0),2);
        for i=1:numel(u_0)
            tmp(i,:) = [u_0(i),v_0];
        end
        u_sol = [tmp(:,1)];
        v_sol = [tmp(:,2)];
    else
        v_sol = [];
        u_sol = [];
    end
else
    v_sol = -r(2,3);
    v_0 = -(r(1,3)+r(2,3)^2)/r(1,2);
    u_sol = v_0;
end
end