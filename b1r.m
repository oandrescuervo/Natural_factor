function f=b1r(xi,eta,r)

x=xi;
y=eta;

f=(y-1+r)*(r-1)+(2*r-1)*x;

x=r; y=0;
fr=(y-1+r)*(r-1)+(2*r-1)*x;
f=f/fr;
