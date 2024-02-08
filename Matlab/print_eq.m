function print_eq(d0, d1, eig_vals, offsets, const)
s = "";
str_d0 = [d0,d1];
if eig_vals(1) ~= 0
    s = strcat(s,sprintf("%+g*(xi%+g)^2",eig_vals(1),-offsets(1)));
else
    s = strcat(s,sprintf("%+g*xi",str_d0(1)));
    str_d0 = str_d0(end);
end
if eig_vals(2) ~= 0
    s = strcat(s,sprintf("%+g*(eta%+g)^2",eig_vals(2),-offsets(2)));
else
    s = strcat(s,sprintf("%+g*eta",str_d0(1)));
    str_d0 = str_d0(end);
end
if eig_vals(3) ~= 0
    s = strcat(s,sprintf("%+g*(zeta%+g)^2",eig_vals(3),-offsets(3)));
else
    s = strcat(s,sprintf("%+g*zeta",str_d0(1)));
end
s = strcat(s,sprintf("%+g",const));
fprintf("%s\n",s)
end