function [v_sol,w_sol] = rank2_3C3_15(r)
    %RANK_3C3_15 Solve the resulting subsystem of equations for the case R15.
    arguments
        r (2,4) {mustBeReal}
    end
    % Obtain root candidates from equations
    w_root = -r(2,4);
    v_root = roots([r(1,1),r(1,3),r(1,2)*r(2,4)^2+r(1,4)-r(2,4)^3]);
    v_sol = v_root(abs(imag(v_root))<10^-10);
    w_sol = repmat(w_root,size(v_sol));
    end