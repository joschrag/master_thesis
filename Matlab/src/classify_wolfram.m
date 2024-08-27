function result = classify_wolfram(c)
arguments
    c (1,10) {mustBeReal}
end
c = c./[1,2,2,1,2,1,2,2,2,1];
e = [c(1),c(2),c(3);...
    c(2),c(4),c(5);...
    c(3),c(5),c(6)];
E = [c(1),c(2),c(3),c(7);...
    c(2),c(4),c(5),c(8);...
    c(3),c(5),c(6),c(9);...
    c(7),c(8),c(9),c(10)];
p_3 = rank(e);
p_4 = rank(E);
d = det(E);
eigs = eig(e);
eigs = (abs(eigs) > 10^-10).*eigs;

nzero_eigs = eigs(eigs~=0);
if all(eigs== 0) 
    error("equation must contain quadratic terms.")
end
k = logical(all(sign(nzero_eigs) == sign(nzero_eigs(1))));

switch p_3
    case 1
        if p_4 == 1
            %one plane
            result = 13;
        elseif p_4 ==  2
            %par planes
            result = 12;
        elseif p_4 ==  3
            %par cylinder
            result = 11;
        else
            result=0;
        end
    case 2
        if p_4 == 2
            if k == 0
                %intersect planes
                result = 14;
            elseif k == 1
                %line
                result = 9;
            else
                result=0;
            end
        elseif p_4 == 3
            if k == 0
                %hyp cylinder
                result = 10;
            elseif k == 1
                %ell cylinder
                result = 8;
            else
            error("");
            end
        elseif p_4 == 4
            if k == 0
                %hyp paraboloid
                result = 7;
            elseif k == 1
                %ell paraboloid
                result = 6;
            else
            result=0;
            end
        else
            result=0;
        end
    case 3
        if p_4 == 3
            if k == 1
                %point solution
                result = 1;
            elseif k == 0
                %ell cone
                result = 2;
            end
        elseif p_4 == 4
            if sign(d) == 1
                if k == 1
                    result=0;
                elseif k == 0
                    %hyperboloid1
                    result = 4;
                else
                    result=0;
                end
            elseif sign(d) == -1
                if k == 1
                    %ellipsoid
                    result = 3;
                elseif k == 0
                    %hyperboloid2
                    result = 5;
                else
                   result=0;
                end
            else
                result=0;
            end
        else
            result=0;
        end
end
end