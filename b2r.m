function f=b2r(xi,eta,r)

x=xi;
y=eta;

f=(y-1+r)*(r)+(1-r)*x;

x=1-r; y=r;
fr=(y-1+r)*(r)+(1-r)*x;
f=f/fr;
