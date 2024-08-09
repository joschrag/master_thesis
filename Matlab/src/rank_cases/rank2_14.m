function [u_sol,v_sol] = rank2_14(r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (2,3) {mustBeReal}
end
u_sol = [];
v_sol = [];
v_0 = -r(2,3);
u_0 = roots([-1,-r(1,2),-r(1,1)*r(2,3)^2-r(1,3)]);
u_0 = u_0(imag(u_0)==0);
if ~isempty(u_0) && ~isempty(v_0)
    tmp = zeros(numel(v_0)*numel(u_0),2);
    for i=1:numel(u_0)
        for j=1:numel(v_0)
            tmp((i-1)*numel(v_0)+j,:) = [u_0(i),v_0(j)];
        end
    end
    u_sol = tmp(:,1);
    v_sol = tmp(:,2);
end
end