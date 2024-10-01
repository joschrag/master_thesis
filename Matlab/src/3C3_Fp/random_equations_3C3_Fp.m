function [c] = random_equations_3C3_Fp(coeff_var,prime)
%RANDOM_EQUATIONS4C3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    coeff_var (1,1) {mustBeMember(coeff_var,[1,2,3])}
    prime (1,1) {mustBePrimeOrZero,mustBePositive};
end
q = [...
    1,7,8,11,15,17,20;...
    1,3,7,13,14,18,20;...
    1,2,10,12,16,19,20;...
    ];
e = [7,8,15;1,3,14;1,2,12];
c = zeros(3,20);
while all(c(:,20) == 0)
    c(:,q(coeff_var,:)) = randi(prime,3,size(q,2))-1;
end
c(:,e(coeff_var,:)) = diag([-1,-1,-1]);
end

