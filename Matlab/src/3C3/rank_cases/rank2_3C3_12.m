function [v_sol,w_sol] = rank2_3C3_12(r,complex)
    %RANK_3C3_12 Solve the resulting subsystem of equations for the case R12.
    arguments
        r (2,4) 
        complex (1,1) {mustBeNumericOrLogical} = false;
    end
    v = sym("v","real");
    w = sym("w","real");
    v_sol = [];
    w_sol = [];
    eq = r(2,1)*w^2+r(2,2)*v+r(2,3)*w+r(2,4)+v^2;
    % Obtain root candidates from equations
    if r(1,2) ~= 0
        v_0 = -(w^3+r(1,1)*w^2+r(1,3)*w+r(1,4))/r(1,2);
        w_root = roots(coeffs(subs(eq,v,v_0),w,"All"));
        if ~complex
        w_root =  w_root(abs(imag(w_root))<10^-10);
        end
        w_sol = w_root;
        v_sol = subs(v_0,w,w_sol);
    else
        w_root = roots([1,r(1,1),r(1,3),r(1,4)]);
        if ~complex
        w_root =  w_root(abs(imag(w_root))<10^-10);
        end
        for sol = w_root'
            v_root = roots(coeffs(subs(eq,w,sol),v,"All"));
            if ~complex
                v_root = v_root(abs(imag(v_root))<10^-10);
            end
            if ~isempty(w_root) & ~isempty(v_root)
                w_sol = [w_sol;repmat(sol,size(v_root))];
                v_sol = [v_sol;v_root];
            end
        end
    end
    end