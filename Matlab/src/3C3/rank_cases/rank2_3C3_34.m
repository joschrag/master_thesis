function [v_sol,w_sol] = rank2_3C3_34(r,complex)
    %RANK_3C3_34 Solve the resulting subsystem of equations for the case R34.
    arguments
        r (2,4) 
        complex (1,1) {mustBeNumericOrLogical} = false;
    end
    w_root = roots([1,r(1,3),r(1,4)]);
    if ~complex
    w_root = w_root(abs(imag(w_root))<10^-10);
    end
    w_sol = w_root;
    v_sol = -r(2,3)*w_sol-r(2,4);
    end