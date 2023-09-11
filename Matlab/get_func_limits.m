function [limits] = get_func_limits(p)
p = p(find(p~=0,1):end);
is_even_degree = logical(mod(numel(p),2)) + 1;
p_sign = sign(p(1));

limits = [(-1)^(is_even_degree)*p_sign,p_sign];

end