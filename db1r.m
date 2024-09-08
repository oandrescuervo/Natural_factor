function [fx,fy]=db1r(r)


fx=2*r-1;
fy=r-1;

%f=(y-1+r)*(r-1)+(2*r-1)*x;

x=r; y=0;
fr=(y-1+r)*(r-1)+(2*r-1)*x;

fx=fx/fr;
fy=fy/fr;