function [v_sol,w_sol] = rank2_3C3_24(r,complex)
    %RANK_3C3_23 Solve the resulting subsystem of equations for the case R23.
    arguments
        r (2,4) 
        complex (1,1) {mustBeNumericOrLogical} = false;
    end
    v = sym("v","real");
    w = sym("w","real");
    eq = r(1,2)*w^2+r(1,3)*w+r(1,4)+v^2;
    v_root = -r(2,3)*w-r(2,4);
    w_root = roots(coeffs(subs(eq,v,v_root),v,"All"));
    if ~complex
        w_root = w_root(abs(imag(w_root))<10^-10);
    end
    w_sol = w_root;
    v_sol = repmat(v_root,size(w_sol));
    end