function result = classify_wolfram(c)
arguments
    c (1,10) {mustBeReal}
end
c = c./[1,1,1,2,2,2,2,2,2,1];
e = [c(1),c(4),c(5);...
    c(4),c(2),c(6);...
    c(5),c(6),c(3)];
E = [c(1),c(4),c(5),c(7);...
    c(4),c(2),c(6),c(8);...
    c(5),c(6),c(3),c(9);...
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
            fprintf("one plane\n")
            result = 13;
        elseif p_4 ==  2
            fprintf("par planes\n")
            result = 12;
        elseif p_4 ==  3
            fprintf("par cylinder\n")
            result = 11;
        else
            result=0;
        end
    case 2
        if p_4 == 2
            if k == 0
                fprintf("intersect planes\n")
                result = 14;
            elseif k == 1
                fprintf("line\n")
                result = 9;
            else
                result=0;
            end
        elseif p_4 == 3
            if k == 0
                fprintf("hyp cylinder\n")
                result = 10;
            elseif k == 1
                fprintf("ell cylinder\n")
                result = 8;
            else
            error("");
            end
        elseif p_4 == 4
            if k == 0
                fprintf("hyp paraboloid\n")
                result = 7;
            elseif k == 1
                fprintf("ell paraboloid\n")
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
                fprintf("point solution\n")
                result = 1;
            elseif k == 0
                fprintf("ell cone\n")
                result = 2;
            end
        elseif p_4 == 4
            if sign(d) == 1
                if k == 1
                    result=0;
                elseif k == 0
                    fprintf("hyperboloid1\n")
                    result = 4;
                else
                    result=0;
                end
            elseif sign(d) == -1
                if k == 1
                    fprintf("ellipsoid\n")
                    result = 3;
                elseif k == 0
                    fprintf("hyperboloid2\n")
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