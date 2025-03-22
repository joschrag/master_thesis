function [v_sol,w_sol] = rank2_3C3_23(r)
    %RANK_3C3_23 Solve the resulting subsystem of equations for the case R23.
    arguments
        r (2,4) {mustBeReal}
    end
    v = sym("v","real");
    w = sym("w","real");
    v_sol = [];
    w_sol = [];
    
    % Obtain root candidates from equations
    if r(1,3) ~= 0
        eq = r(2,2)*v+r(2,3)*w+r(2,4)+w^2;
        w_0 = -(v^2+r(1,2)*v+r(1,4))/r(1,3);
        v_root = roots(coeffs(subs(eq,w,w_0),w,"All"));
        v_sol =  v_root(abs(imag(v_root))<10^-10);
        w_sol = subs(w_0,v,v_sol);
    elseif r(2,2) ~= 0
        eq = r(1,2)*v+r(1,3)*w+r(1,4)+v^2;
        v_0 = -(w^2+r(2,3)*w+r(2,4))/r(2,2);
        w_root = roots(coeffs(subs(eq,v,v_0),w,"All"));
        w_sol =  w_root(abs(imag(w_root))<10^-10);
        v_sol = subs(v_0,w,w_sol);
    else
        v_root = roots([1,r(1,2),r(1,4)]);
        w_root = roots([1,r(2,3),r(2,4)]);
        v_root =  v_root(abs(imag(v_root))<10^-10);
        w_root =  w_root(abs(imag(w_root))<10^-10);
        if ~isempty(v_root) & ~isempty(w_root)
            [v_sol,w_sol] = meshgrid(v_root,w_root);
            v_sol = v_sol(:);
            w_sol = w_sol(:);
        end
    end
    end