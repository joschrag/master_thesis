function [v_sol,w_sol] = rank2_5C3_24_Fp(r,prime)
%RANK2_5C3_24_FP Solve the resulting subsystem of equations for the case R24.
arguments
    r (2,3) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePrimeOrZero,mustBePositive}
end
% Obtain solutions from equations based on coefficients
if r(1,2) == 0
    if FF(r(1,3),prime) == -FF(r(2,3),prime)^2
        w_sol = 1:prime;
        v_sol = FF(repmat(-r(2,3),size(w_sol)),prime).value;
    else
        w_sol = [];
        v_sol = [];
    end
else
    w_sol = FF(-r(2,3),prime).value;
    v_sol = -FF(r(1,3)+r(2,3)^2,prime)*FF(r(1,2),prime)^(-1);
    v_sol = v_sol.value;
end
end