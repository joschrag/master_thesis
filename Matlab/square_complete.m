function [offset,const] = square_complete(a,b,c)

offset = (b/(2*a));
const = c-(b^2/(4*a));

end