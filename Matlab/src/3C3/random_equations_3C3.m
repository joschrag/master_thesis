function [c] = random_equations_3C3(coeff_var)
%RANDOM_EQUATIONS4C3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    coeff_var (1,1) {mustBeMember(coeff_var,[1,2,3])}
end
q = [...
    1,7,8,11,15,17,20;...
    1,3,7,13,14,18,20;...
    1,2,10,12,16,19,20;...
    ];
c = zeros(3,20);
while all(c(:,20) == 0)
    c(:,q(coeff_var,:)) = randi(100,3,size(q,2))-50;
end
end

