function [v_sol,w_sol] = rank2_5C3_24(r)
%RANK2_5C3_24 Solve the resulting subsystem of equations for the case R24.
arguments
    r (2,3) {mustBeReal}
end
% Obtain solutions from equations based on coefficients
if r(1,2) == 0
    if r(1,3) == -r(2,3)^2
        w_0 = -r(2,3);
        v_0 = -10:0.01:10;
        [v_sol,w_sol] = meshgrid(v_0,w_0);
    else
        w_sol = [];
        v_sol = [];
    end
else
    w_sol = -r(2,3);
    v_0 = -(r(1,3)+r(2,3)^2)/r(1,2);
    v_sol = v_0;
end
end