function [fx,fy]=db3r(r)


%f=(y-r)*(1-2*r)-(r)*(x-1+r);
fx=-r;
fy=1-2*r;
x=0; y=1-r;
fr=(y-r)*(1-2*r)-(r)*(x-1+r);
fx=fx/fr;
fy=fy/fr;