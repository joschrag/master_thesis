function [limits] = get_func_limits(p)
arguments
    p (:,:) {mustBeReal};
end
[n,m] = size(p);
c = zeros(n,1);
s = zeros(n,1);
for i=1:n
c(i) = find(p(i,:)~=0,1,"first");
s(i) = sign(p(i,c(i)));
end
is_even_degree = logical(mod(m-c,2));
limits = [(-1).^(is_even_degree).*s,s]';
end