function kappa=coefficient_contraste(x1,x2,va)
kappa=zeros(size(x1));

for i=1:size(x1,1)
  x=x1(i);
  y=x2(i);
  
  if (0.1<x) && (x<0.2) && (0.1<y) && (y<0.9)
      
      kappa(i)=exp(va(1));
      
  elseif (0.3<x) && (x<0.4) && (0.1<y) && (y<0.9)
      
      kappa(i)=exp(va(2));
      
  elseif (0.5<x) && (x<0.6) && (0.1<y) && (y<0.9)
      
      kappa(i)=exp(va(3));
      
 elseif (0.7<x) && (x<0.8) && (0.1<y) && (y<0.9)
      
      kappa(i)=exp(va(4));
           
  else
      
      kappa(i)=1;
      
  end
  
end
      
end



%kappa=1+10000*(2+sin(8*pi*(x1.^2+x2)));
%kappa=1+10*(2+sin(8*pi*(x1.^2+x2)));


%kappa=exp(randn(mesh.ne,1));

