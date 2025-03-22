function [v_sol,w_sol] = rank2_3C3_25(r,complex)
    %RANK_3C3_23 Solve the resulting subsystem of equations for the case R23.
    arguments
        r (2,4) 
        complex (1,1) {mustBeNumericOrLogical} = false;
    end
    w_root = -r(2,4);
    if r(1,3) ~= 0
        v_root = roots([1,r(1,3),r(1,2)*r(2,4)^2+r(1,4)]);
        if ~complex
            v_root = v_root(abs(imag(v_root))<10^-10);
        end
        v_sol = v_root;
        w_sol = repmat(w_root,size(v_sol));
    else
        v_root = roots([1,0,r(1,2)*r(2,4)^2+r(1,4)]);
        if ~complex
        v_root = v_root(abs(imag(v_root))<10^-10);
        end
        v_sol = v_root;
        w_sol = repmat(w_root,size(v_sol));
    end
    end