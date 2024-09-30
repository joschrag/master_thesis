function [u_sol,v_sol] = rank2_5C3_24_Fp(r,prime)
%RANK2_5C3_24_FP Solve the resulting subsystem of equations for the case R24.
arguments
    r (2,3) {mustBeReal,mustBeInteger}
    prime (1,1) {mustBePrime}
end
% Obtain solutions from equations based on coefficients
if r(1,2) == 0
    if FF(r(1,3),prime) == -FF(r(2,3),prime)^2
        v_sol = 1:prime;
        u_sol = FF(repmat(-r(2,3),size(v_sol)),prime).value;
    else
        v_sol = [];
        u_sol = [];
    end
else
    v_sol = FF(-r(2,3),prime).value;
    u_sol = -FF(r(1,3)+r(2,3)^2,prime)*FF(r(1,2),prime)^(-1);
    u_sol = u_sol.value;
end
end