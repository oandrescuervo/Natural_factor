clear all
clc

ax=0;bx=1;ay=0;by=1;

n=20;  %number of elements in each direction
R=1000;  %Número de realizaciones MC
K=30;  %Número de términos en la serie KL


intx=[ax,bx];inty=[ay,by];

%mesh
[M,B,v,mesh,h,free,MC]=createmeshr(n,intx,inty);
%v guarda los vértices (primera columna coordenada en x y segunda columna coordenada en y)
%M guarda las coordenada de los vértices de cada elemento
%MC Puntos medios de cada elemento



nvel=size(v,1);  %Número de vértices
[E,DOF,mesh,efree]=edgesenum(M,B,v,mesh,h,free,MC,n);   %%%función para numerar los lados de los elementos - está un poco desorganizado




% Global Neumann matrix for Crouzeixâ€?Raviart (CR) element
r=0.5; % CR mid point based.
%[AnCR,bCR]=Nmatrixr(M,DOF,v,mesh,kappa,r);
bCR=rhsCR(M,DOF,v,mesh,r);                       %%%Cálculo del bCR

%Global Neumann matrix of Lagrange (L) elements
%[AnL,bL]=NmatrixrP1(M,v,mesh,kappa);
bL=rhsP1(M,v,mesh);                              %%%Cálculo del bL



%%%La matriz de elementos finitos se monta con el gradiente de las
%%%funciones base en el centro del triángulo %%%% Prestar atención a esto
%%%para CR

%Matrix G_L for partial_x partial_y
[GG1,GG2]=GmatrixP1(M,v,mesh);
list2=1:mesh.nv-1; % we need to take out 1 degree of freedom!
GL=[GG2;-GG1];              %%%Está en la forma rotacional para generar independencia
clear GG2 GG1
%Matrix G_CR for partial_x partial_y
[G1,G2]=Gmatrixr(M,v,mesh,DOF,r);
GCR=[G1;G2];
clear G1 G2
% epsilon for test
myeps=1;


%%%efree es el listado de las aristas interiores
%%%list2 son todos los vértices menos el último
%%%free es el listado de los vértices interiores (no se usa en esta aplicación)
Gepsilon=[GCR(:,efree),myeps*GL(:,list2)];
G0=[GCR(:,efree),GL(:,list2)];                  %%% Matriz general G (H en artículo)
G1=[GCR(:,efree),-GL(:,list2)];                  %%% Matriz general G (H en artículo)

% Inversa de matriz H - precomputo
G0_inv=inv(G0);
G1_inv=inv(G1);



bepsilon=[bCR(efree);myeps*bL(list2)];
b=[bCR(efree);bL(list2)];                      %%% Vector del lado derecho b



%%%%------------------------------------
%%%% DEPENDECE ON KAPPA
%%%%


% Parámetros función de covarianza (solo si se define en VP usar la Matern)
nu=0.5;
l=1;


% Función que calcula los valores propios
[V,e]=VP(MC,nu,l);
% V: vectores propios
% e: valores propios 

% Números aleatorios
randn('seed',1);
va=randn(K,R); 
    
for j=1:R
   
% Coeficiente c expansión KL     
c=0;
for i=1:K
    c=c+sqrt(e(i))*va(i,j)*V(:,i);
end

% coefficient
kappa=exp(c);


% Diagonal with the cofficient in partial_x and partial_y
D=sparse(2*mesh.ne,2*mesh.ne);
D_inv=sparse(2*mesh.ne,2*mesh.ne);
for iaux=1:mesh.ne
    D(iaux,iaux)=kappa(iaux);
    D(mesh.ne+iaux,mesh.ne+iaux)=kappa(iaux);
    D_inv(iaux,iaux)=1/kappa(iaux);
    D_inv(mesh.ne+iaux,mesh.ne+iaux)=1/kappa(iaux);
    coef(iaux)=kappa(iaux);
end

contraste(j)=max(coef)/min(coef);



% Block diagonal matrix 
%zaux=sparse(length(list2),mesh.nedgesf);
%AL=GL(:,list2)'*D*GL(:,list2);
%ACR=GCR(:,efree)'*D*GCR(:,efree);
%ABDepsilon=[ACR,zaux';zaux,myeps*AL];

%ABD=[ACR,zaux';zaux,AL];                  %%%Matriz diagonal por bloques A^tilde=[ACR,0;0,AL]

% Matrices nuevo sistema
A0=sparse(2*mesh.ne,2*mesh.ne);
A1=sparse(2*mesh.ne,2*mesh.ne);
M0=sparse(2*mesh.ne,2*mesh.ne);
M1=sparse(2*mesh.ne,2*mesh.ne);


A0=G0'*D*G0;
A1=G1'*D*G1;
M0=G0_inv*D_inv*G0_inv';
M1=G1_inv*D_inv*G1_inv';
M_05=0.5*M0+0.5*M1;

A_sys=A1*M_05*A0;

% Multiplication of 3 square matrices
MatrizH1=GCR(:,efree)'*GCR(:,efree);
%MatrixM=G'*D*G;                                %%%Matriz M
%MatrixM2=Gepsilon'*D*G;
%MatrixM3=G'*D*Gepsilon;

% free degrees of freedom for block diagonal system
tolfree=[efree,mesh.nedges+(1:length(list2))];

u1=b*0;u2=b*0;u3=b*0;u4=b*0;u5=b*0;


% Precond. using LU ---- Solución del método Iterativo
% Método anterior
[xLU, errorLU, iterLU, flagLU, lambdamaxLU, condnumberLU] =ApcgCR(b*0, b,b, 20000, 0.000001, A_sys,A1);


NCond(j)=condnumberLU;
Iter(j)=iterLU;




end




media_NCond=mean(NCond)
varianza_NCond=var(NCond)
media_Iter=mean(Iter)
varianza_Iter=var(Iter)
media_contraste=mean(contraste)
varianza_contraste=var(contraste)


%Error=sqrt((solCR2-sol_prom)'*MatrizH1*(solCR2-sol_prom));




%Error

%%GRÁFICA VARIANDO L

%CondNumber(cont,2)=mean(NCond(:,1));
%CondNumber(cont,1)=l;
%plot(CondNumber(:,1),CondNumber(:,2))


% hist(NCond)
% title('Matrix condition histogram')
% xlabel('Matrix condition number')
% ylabel('frequency')

