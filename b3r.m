function f=b3r(xi,eta,r)

x=xi;
y=eta;

f=(y-r)*(1-2*r)-(r)*(x-1+r);

x=0; y=1-r;
fr=(y-r)*(1-2*r)-(r)*(x-1+r);
f=f/fr;
