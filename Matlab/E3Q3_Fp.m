function [result] = E3Q3_Fp(c,prime)
%E3Q3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    c (3,10) {mustBeInteger} = [1, 2, 1, 0, 1, 1, 2, 0, 0, -4;...
        1, 1, 3, 2, 1, -1, 2, 0, 2, -5;...
        1, -1, 1, 1, 0, 1, -1, 0, 3, -6;
        ];
    prime (1,1) {mustBeInteger} = nextprime(6);
end
x = sym("x","integer");
y = sym("y","integer");
z = sym("z","integer");
v = [x^2,y^2,z^2,x*y,x*z,y*z,x,y,z,1]';
A = -FF([c(:,2),c(:,3),c(:,6)],prime);
B = -FF([c(:,1),c(:,3),c(:,5)],prime);
C = -FF([c(:,1),c(:,2),c(:,4)],prime);
[m,~] = size(c);
r_A = rank(A);
r_B = rank(B);
r_C = rank(C);
if any([r_A,r_B,r_C] == 3)
    if r_A == 3
        fprintf("Using P(x)\n")
        Q = A;
        P = FF([c(:,4).*x+c(:,8),c(:,5).*x+c(:,9),c(:,1).*x^2+c(:,7).*x+c(:,10)],prime);
        lin_vars = FF([y;z;1],prime);
        p_var = x;
    elseif r_B == 3
        fprintf("Using P(y)\n")
        Q = B;
        P = FF([c(:,4).*y+c(:,7),c(:,6).*y+c(:,9),c(:,2).*y^2+c(:,8).*y+c(:,10)],prime);
        lin_vars = FF([x;z;1],prime);
        p_var = y;
    elseif r_C == 3
        fprintf("Using P(z)\n")
        Q = C;
        P = FF([c(:,5).*z+c(:,7),c(:,6).*z+c(:,8),c(:,3).*z^2+c(:,9).*z+c(:,10)],prime);
        lin_vars = FF([x;y;1],prime);
        p_var = z;
    end
    P2 = inv(Q)*P;

    quad_1 = FF(P2.value(1,:),prime)*lin_vars;
    quad_2 = FF(P2.value(2,:),prime)*lin_vars;
    mixed = FF(P2.value(3,:),prime)*lin_vars;

    identities = [(quad_1)*FF(lin_vars.value(2),prime) == (mixed)*FF(lin_vars.value(1),prime);...
        (mixed)*FF(lin_vars.value(2),prime) == (quad_2)*FF(lin_vars.value(1),prime);...
        (mixed)*(mixed) == quad_1*quad_2];
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
    eq = simplify(lhs(eq)-rhs(eq));
    [A_,b] = equationsToMatrix(eq,lin_vars.value(1:2));
    A2 = FF([A_,-b],prime);
    pol = coeffs(det(A2),p_var,"All");
    [~,x_sol,~] = gfroots(fliplr(pol),1,prime);
    result = [];
    for x_0 = x_sol
        result(1) = x_0;
        M = subs(A2,x,x_0);
        b_vec = -M;
        [result(2:3),~] = gflineq(double(M.value(1:2,1:2)),double(b_vec.value(1:2,3)),prime);
    end

    for i=1:size(result,1)
        n = norm(mod(eval(subs(c*v,[x,y,z],result(i,:))),prime));
        eval_str = join(repmat("\n%g",[1,m]),"");
        line_break = join(repmat("=",[1,10]),"");
        res_str = join(repmat("%i",[1,m]),",");
        res_str = join(["(",res_str,")"]);
        result_str = sprintf("Lösung %%i %s-> \nNorm: %%i%s\n%%s\n",res_str,eval_str);
        fprintf(result_str,i,result(i,:),n,mod(eval(subs(c*v,[x,y,z],result(i,:))),prime),line_break);
    end


end
for i =0:prime-1
    for k =0:prime-1

        for j =0:prime-1
            if  ~any(mod(eval(subs(c*v,[x,y,z],[i,k,j])),prime))
            fprintf("%i %i %i | %i %i %i\n",i,k,j,mod(eval(subs(c*v,[x,y,z],[i,k,j])),prime))
            
            end
        end
    end
end
end

