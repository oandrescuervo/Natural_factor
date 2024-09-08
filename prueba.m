
ax=0;bx=1;ay=0;by=1;

n=3; %number of elements in each direction
intx=[ax,bx];inty=[ay,by];


[M,B,v,mesh,h,free,MC]=createmeshr(n,intx,inty);

xi=[ power(3,-1), (6+sqrt(15))/21 ,(9-2*sqrt(15))/21,(6+sqrt(15))/21,(6-sqrt(15))/21 ,(9+2*sqrt(15))/21,(6-sqrt(15))/21];
eta=[power(3,-1),(6+sqrt(15))/21,(6+sqrt(15))/21,(9-2*sqrt(15))/21,(6- sqrt(15))/21,(6-sqrt(15))/21, (9+2*sqrt(15))/21];


for i=1:mesh.ne
T=M(i,:);
for j=1:3
    x(j)=v(T(j),1);
    y(j)=v(T(j),2);
end

x1=x(1); x2=x(2); x3=x(3);
y1=y(1); y2=y(2); y3=y(3);

txi=x(1)+xi*(x(2)-x(1))+eta*(x(3)-x(1));
teta=y(1)+xi*(y(2)-y(1))+eta*(y(3)-y(1));


%function_f1=f(txi,teta);


%x=xi;
%y=eta;

f=(y-1+r)*(r-1)+(2*r-1)*x;

%x=r; y=0;
%fr=(y-1+r)*(r-1)+(2*r-1)*x;
%f=f/fr;

end

