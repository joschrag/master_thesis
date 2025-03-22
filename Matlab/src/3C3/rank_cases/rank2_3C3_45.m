function [v_sol,w_sol] = rank2_3C3_45(r)
    %RANK_3C3_15 Solve the resulting subsystem of equations for the case R15.
    arguments
        r (2,4) {mustBeReal}
    end
    % Obtain root candidates from equations
    w_sol = -r(2,4);
    v_sol = -r(1,4);
    end