function exact=analytic(x,y)
exact=sin(2*pi*x).*sin(pi*y).*exp(-(x+y));
%exact=sin(2*pi.*x).*sin(pi.*y).*exp(x+y);