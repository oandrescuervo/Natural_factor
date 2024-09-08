function L=localb(T,v,par,r)


for j=1:3
    x(j)=v(T(j),1);
    y(j)=v(T(j),2);
end




omega=[9*power(80,-1), (155+sqrt(15))/2400, (155+sqrt(15))/2400,(155+sqrt(15))/2400,(155-sqrt(15))/2400, (155-sqrt(15))/2400,(155-sqrt(15))/2400];

xi=[ power(3,-1), (6+sqrt(15))/21 ,(9-2*sqrt(15))/21,(6+sqrt(15))/21,(6-sqrt(15))/21 ,(9+2*sqrt(15))/21,(6-sqrt(15))/21];

eta=[power(3,-1),(6+sqrt(15))/21,(6+sqrt(15))/21,(9-2*sqrt(15))/21,(6- sqrt(15))/21,(6-sqrt(15))/21, (9+2*sqrt(15))/21];

den=abs(x(2)*y(1)-x(3)*y(1)-x(1)*y(2)+x(3)*y(2)+x(1)*y(3)-x(2)*y(3));


txi=x(1)+xi*(x(2)-x(1))+eta*(x(3)-x(1));
teta=y(1)+xi*(y(2)-y(1))+eta*(y(3)-y(1));


x1=x(1); x2=x(2);x3=x(3);
y1=y(1); y2=y(2);y3=y(3);
B=[x2-x1, x3-x1;y2-y1,y3-y1];
InvB=inv(B);
ub_x(23)=2*InvB(1,1);ub_y(2)=-2*InvB(1,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ub_x(2)=2*InvB(1,1)+2*InvB(2,1);ub_y(3)=2*InvB(2,1)+2*InvB(2,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ub_x(1)=-2*InvB(2,1);ub_y(1)=-2*InvB(2,2);

den=abs(x(2)*y(1)-x(3)*y(1)-x(1)*y(2)+x(3)*y(2)+x(1)*y(3)-x(2)*y(3));



x1=x(1); x2=x(2);x3=x(3);
y1=y(1); y2=y(2);y3=y(3);

pb(1,:)=b1r(xi,eta,r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pb(2,:)=b2r(xi,eta,r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pb(3,:)=b3r(xi,eta,r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%meanp=0;
%for i=1:3
%    meanp=meanp+u(2*par.nvel1+index(i))*pb(i,:);
%end
function_f1=f(txi,teta);
for i=1:3
   L(i)=sum( omega.*pb(i,:).*function_f1);
%   L=L+sum( transpose(omega).*function_f2);
end

B=[x2-x1, x3-x1;y2-y1,y3-y1];

L=L*det(B);
%den=abs(x(2)*y(1)-x(3)*y(1)-x(1)*y(2)+x(3)*y(2)+x(1)*y(3)-x(2)*y(3));


