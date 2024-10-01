function print_solutions(result,eqs,x,y,z,prime)
arguments
    result (:,3);
    eqs (:,:);
    x (1,1) = sym("x","real");
    y (1,1) = sym("y","real");
    z (1,1) = sym("z","real");
    prime (1,1) {mustBePrimeOrZero} = 0;
end
n = zeros(size(result,1),1);
m = size(eqs,1);
for i=1:size(result,1)
    n(i) = norm(mod(eval(subs(eqs,[x,y,z],result(i,:))),prime));
    eval_str = join(repmat("\n%g",[1,m]),"");
    line_break = join(repmat("=",[1,10]),"");
    res_str = join(repmat("%g",[1,3]),",");
    res_str = join(["(",res_str,")"]);
    result_str = sprintf("Solution %%i %s-> \nNorm: %%i%s\n%%s\n",res_str,eval_str);
    fprintf(result_str,i,result(i,:),n(i),mod(eval(subs(eqs,[x,y,z],result(i,:))),prime),line_break);
end
end