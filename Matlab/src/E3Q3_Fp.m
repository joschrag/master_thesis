function [result] = E3Q3_Fp(c,prime)
%E3Q3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    c (:,10) {mustBeInteger} = [1, 2, 1, 0, 1, 1, 2, 0, 0, -4;...
        1, 1, 3, 2, 1, -1, 2, 0, 2, -5;...
        1, -1, 1, 1, 0, 1, -1, 0, 3, -6;
        ];
    prime (1,1) {mustBeInteger} = nextprime(6);
end
x = sym("x","integer");
y = sym("y","integer");
z = sym("z","integer");
vars = [x,y,z];
v = [x^2,x*y,x*z,y^2,y*z,z^2,x,y,z,1]';
A = -FF([c(:,4),c(:,5),c(:,6)],prime);
B = -FF([c(:,1),c(:,3),c(:,5)],prime);
C = -FF([c(:,1),c(:,2),c(:,4)],prime);
r_A = rank(A);
r_B = rank(B);
r_C = rank(C);
if any([r_A,r_B,r_C] == 3)
    if r_A == 3
        fprintf("Using P(x)\n")
        Q = A;
        P = FF([c(:,2).*x+c(:,8),c(:,3).*x+c(:,9),c(:,1).*x^2+c(:,7).*x+c(:,10)],prime);
        lin_vars = FF([y;z;1],prime);
        p_var = x;
    elseif r_B == 3
        fprintf("Using P(y)\n")
        Q = B;
        P = FF([c(:,2).*y+c(:,7),c(:,5).*y+c(:,9),c(:,4).*y^2+c(:,8).*y+c(:,10)],prime);
        lin_vars = FF([x;z;1],prime);
        p_var = y;
    elseif r_C == 3
        fprintf("Using P(z)\n")
        Q = C;
        P = FF([c(:,3).*z+c(:,7),c(:,5).*z+c(:,8),c(:,6).*z^2+c(:,9).*z+c(:,10)],prime);
        lin_vars = FF([x;y;1],prime);
        p_var = z;
    end

    P2 = l_inv(Q)*P;
    quad_1 = FF(P2.value(1,:),prime)*lin_vars;
    quad_2 = FF(P2.value(3,:),prime)*lin_vars;
    mixed = FF(P2.value(2,:),prime)*lin_vars;

    identities = [(quad_1)*FF(lin_vars.value(2),prime) - (mixed)*FF(lin_vars.value(1),prime);...
        (mixed)*FF(lin_vars.value(2),prime) - (quad_2)*FF(lin_vars.value(1),prime);...
        (mixed)*(mixed) - quad_1*quad_2];
    sub1 = FF.empty(3,0);
    sub2 = FF.empty(3,0);
    sub3 = FF.empty(3,0);

    for i = 1:3
        sub1(i) = collect(identities(i));
    end
    old_vars = [lin_vars.value(1)^2,...
        lin_vars.value(2)^2,...
        lin_vars.value(1)*lin_vars.value(2)];
    new_vars = [quad_1.value,quad_2.value,mixed.value];
    for i = 1:3
        sub2(i) = subs(sub1(i),old_vars,new_vars);
    end
    for i = 1:3
        ex = expand(sub2(i));
        sub3(i) = subs(ex,old_vars,new_vars);
    end
    eq = [sub3(1).value,sub3(2).value,sub3(3).value];
    eq = simplify(eq);
    [A_,b] = equationsToMatrix(eq,lin_vars.value(1:2));
    A2 = FF([A_,-b],prime);
    pol = coeffs(det(A2),p_var,"All");
    p_root = get_gf_root(pol,prime);
    if numel(p_root) == 0
        fprintf("No solutions to the equations.\n")
        result = [];
        return;
    end
    idx = [find(vars==p_var),find(vars==lin_vars.value(1)),find(vars==lin_vars.value(2))];
    result = [];
    for i = 1:numel(p_root)
        p_0 = p_root(i);
        M = subs(A2,p_var,p_0);
        if M.value == zeros(3)
            fprintf("Matrix is singular.\n");
            continue;
        end
        disp(M.value)
        rM = rref(M);
        m = zeros(1,rank(rM));
        for j=1:rank(rM)
            m(j) = find(rM(j,:),1,"first");
        end
        r = vpa(rM(1:rank(rM),setdiff(1:3,m)));
        switch join(string(m),"")
            case "1"
                r_root = (0:prime-1)';
                q_root = FF(-r(1).*r_root-r(2),prime).value;
            case "2"
                q_root = (0:prime-1)';
                r_root = FF(-r(2),prime).value;
            case "12"
                o_sols = FF(-r',prime).value;
                q_root = o_sols(1);
                r_root = o_sols(2);
            otherwise
                fprintf("Matrix equations are inconsistent!");
                continue;
        end
        if ~isempty(q_root)
            for root=[q_root,r_root]'
                result=[result;[p_0,root']];
                result(end,:) = result(end,idx);
            end
        end
    end
    equations = c*v;
    if numel(result) > 0
        print_solutions(result,equations,x,y,z,prime)
    end
end
end

