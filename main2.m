clear all

ax=0;bx=1;ay=0;by=1;

n=40; %number of elements in each direction
R=1000;  %N�mero de realizaciones MC
K=15;    %N�mero de t�rminos en la serie KL

intx=[ax,bx];inty=[ay,by];

%mesh
[M,B,v,mesh,h,free,MC]=createmeshr(n,intx,inty);
%v guarda los v�rtices (primera columna coordenada en x y segunda columna coordenada en y)
%M guarda las coordenada de los v�rtices de cada elemento
%MC Puntos medios de cada elemento



nvel=size(v,1);  %N�mero de v�rtices
[E,DOF,mesh,efree]=edgesenum(M,B,v,mesh,h,free,MC,n);   %%%funci�n para numerar los lados de los elementos - est� un poco desorganizado




% Global Neumann matrix for Crouzeix�?Raviart (CR) element
r=0.5; % CR mid point based.
%[AnCR,bCR]=Nmatrixr(M,DOF,v,mesh,kappa,r);
bCR=rhsCR(M,DOF,v,mesh,r);                       %%%C�lculo del bCR

%Global Neumann matrix of Lagrange (L) elements
%[AnL,bL]=NmatrixrP1(M,v,mesh,kappa);
bL=rhsP1(M,v,mesh);                              %%%C�lculo del bL



%%%La matriz de elementos finitos se monta con el gradiente de las
%%%funciones base en el centro del tri�ngulo %%%% Prestar atenci�n a esto
%%%para CR

%Matrix G_L for partial_x partial_y
[GG1,GG2]=GmatrixP1(M,v,mesh);
list2=1:mesh.nv-1; % we need to take out 1 degree of freedom!
GL=[GG2;-GG1];              %%%Est� en la forma rotacional para generar independencia
clear GG2 GG1
%Matrix G_CR for partial_x partial_y
[G1,G2]=Gmatrixr(M,v,mesh,DOF,r);
GCR=[G1;G2];
clear G1 G2
% epsilon for test
myeps=1;

%%%efree es el listado de las aristas interiores
%%%list2 son todos los v�rtices menos el �ltimo
%%%free es el listado de los v�rtices interiores (no se usa en esta aplicaci�n)
Gepsilon=[GCR(:,efree),myeps*GL(:,list2)];
G=[GCR(:,efree),GL(:,list2)];                  %%% Matriz general G

bepsilon=[bCR(efree);myeps*bL(list2)];
b=[bCR(efree);bL(list2)];                      %%% Vector del lado derecho b


%%%Se puede poner lo que quiera en bL por ser el problema auxiliar

%LU factorization of G
[LG,UG]=lu(G);


%%%%------------------------------------
%%%% DEPENDECE ON KAPPA
%%%%


%%%Cambiar esta parte implementando las funciones de covarianza





%PAR�METROS FUNCI�N COVARIANZA
p=0.5;
l=1;
[V,VPDiag,e]=VP(MC,size(MC,1),p,l);


%FOR PARA REALIZAR LOS ERRORES A MEDIDA QUE CRECE K

 
    
    randn('seed',1);
    va=randn(K,R); 

    
    sol=zeros(3*n^2-2*n,R);
    
for j=1:R

c=0;
for i=1:K
    c=c+sqrt(e(i))*va(i,j)*V(:,size(MC,1)+1-i);
end



%coefficient

%kappa=coefficient(MC(:,1),MC(:,2));
kappa=exp(c);


% Diagonal with the cofficient in partial_x and partial_y
D=sparse(2*mesh.ne,2*mesh.ne);
for iaux=1:mesh.ne
    D(iaux,iaux)=kappa(iaux);
    D(mesh.ne+iaux,mesh.ne+iaux)=kappa(iaux);
    coef(iaux)=kappa(iaux);
end

contraste(j)=max(coef)/min(coef);


end



mean(contraste)
var(contraste)

%Error=sqrt((solCR2-sol_prom)'*MatrizH1*(solCR2-sol_prom));




%Error

%%GR�FICA VARIANDO L

%CondNumber(cont,2)=mean(NCond(:,1));
%CondNumber(cont,1)=l;
%plot(CondNumber(:,1),CondNumber(:,2))


% hist(NCond)
% title('Matrix condition histogram')
% xlabel('Matrix condition number')
% ylabel('frequency')

