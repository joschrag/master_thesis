function [roots] = get_gf_root(poly,order)
arguments
    poly (1,:)
    order (1,1) {mustBeInteger,mustBeNonnegative}
end
if isa(poly,"sym")
    poly = double(poly);
end
if ~isprime(order)
    error("FF:order","Order must be prime.")
end
if isempty(poly)
    roots = 0:order-1;
else
prim_pol = gfprimfd(1,"min",order);
q = fliplr(mod(poly,order));
[~,roots] = gfroots(q,prim_pol,order);
end
end
