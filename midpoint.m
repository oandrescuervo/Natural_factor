function [valx,valy]=midpoint(v,E,efree)

nv=size(v,1);
vertx=sparse(nv,nv);
verty=sparse(nv,nv);
for i=1:nv
    for j=1:nv
    vertx(i,j)=(v(i,1)+v(j,1))/2;  
    verty(i,j)=(v(i,2)+v(j,2))/2;    
    end
end

[row,col]=find(E);

Ex=sparse(nv,nv);
Ey=sparse(nv,nv);

for i=1:size(row)
  Ex(row(i),col(i))=vertx(row(i),col(i)); 
  Ey(row(i),col(i))=verty(row(i),col(i)); 
end

nvi=size(efree,2);
valx=zeros(nvi,1);
valy=zeros(nvi,1);

%bp = waitbar(0,'Calculando')

for i=1:nvi
    
  % xx=i/M;
  % waitbar(xx,bp)
    
   [row2,col2]=find(E==efree(i));
   valx(i)=Ex(row2(1),col2(1));
   valy(i)=Ey(row2(1),col2(1));
end





    