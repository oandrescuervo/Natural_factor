function [fx,fy]=db2r(r)

%f=(y-1+r)*(r)+(1-r)*x;
fx=1-r;
fy=r;
x=1-r; y=r;
fr=(y-1+r)*(r)+(1-r)*x;
fx=fx/fr;
fy=fy/fr;
