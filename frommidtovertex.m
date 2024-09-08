function vv=frommidtovertex(T,vm,v)

x=zeros(3,1);
y=zeros(3,1);
L=zeros(3,1);


for j=1:3
    x(j)=v(T(j),1);
    y(j)=v(T(j),2);
end





xi=[0,1,0];

eta=[0,0,1];




pb(1,:)=b1r(xi,eta,.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pb(2,:)=b2r(xi,eta,.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pb(3,:)=b3r(xi,eta,.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%meanp=0;
%for i=1:3
%    meanp=meanp+u(2*par.nvel1+index(i))*pb(i,:);
%end

vv=vm*0;
for i=1:3
    vv=vv+vm(i)*pb(i,:)';
end



