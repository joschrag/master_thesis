function [v_sol,w_sol] = rank2_3C3_13(r)
    %RANK_3C3_12 Solve the resulting subsystem of equations for the case R12.
    arguments
        r (2,4) {mustBeReal}
    end
    v = sym("v","real");
    w = sym("w","real");
    v_sol = [];
    w_sol = [];
    eq = r(1,1)*v^2+r(1,2)*v+r(1,3)*w+r(1,4)+w^3;
    % Obtain root candidates from equations
    if r(2,2) ~= 0
        v_0 = -(w^3+w^2+r(1,3)*w+r(1,4))/r(2,2);
        w_root = roots(coeffs(subs(eq,v,v_0),v,"All"));
        w_sol =  w_root(abs(imag(w_root))<10^-10);
        v_sol = subs(v_0,w,w_sol);
    else
        w_root = roots([1,r(2,3),r(2,4)]);
        w_root =  w_root(abs(imag(w_root))<10^-10);
        for sol = w_root'
            v_root = roots(coeffs(subs(eq,w,sol),v,"All"));
            if ~isempty(w_root(abs(imag(w_root))<10^-10)) & ~isempty(v_root(abs(imag(v_root))<10^-10))
                tmp_v = v_root(abs(imag(v_root))<10^-10);
                w_sol = [w_sol;repmat(sol,size(tmp_v))];
                v_sol = [v_sol;tmp_v];
            end
        end
    end
    end