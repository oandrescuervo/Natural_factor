clear all
clc

ax=0;bx=1;ay=0;by=1;

%n=10;  %number of elements in each direction
R=100;  %Número de realizaciones MC
K=10;  %Número de términos en la serie KL

intx=[ax,bx];inty=[ay,by];

for w=1:3

  n=15+5*w;

%mesh
[M,B,v,mesh,h,free,MC]=createmeshr(n,intx,inty);
%v guarda los vértices (primera columna coordenada en x y segunda columna coordenada en y)
%M guarda las coordenada de los vértices de cada elemento
%MC Puntos medios de cada elemento



nvel=size(v,1);  %Número de vértices
[E,DOF,mesh,efree]=edgesenum(M,B,v,mesh,h,free,MC,n);   %%%función para numerar los lados de los elementos - está un poco desorganizado




% Global Neumann matrix for Crouzeixâ??Raviart (CR) element
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


% Parámetros función de covarianza
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
end


% Block diagonal matrix 
zaux=sparse(length(list2),mesh.nedgesf);
AL=GL(:,list2)'*D*GL(:,list2);
ACR=GCR(:,efree)'*D*GCR(:,efree);
%ABDepsilon=[ACR,zaux';zaux,myeps*AL];

ABD=[ACR,zaux';zaux,AL];                  %%%Matriz diagonal por bloques A^tilde=[ACR,0;0,AL]

% Matrices aproximación 1
A0=sparse(2*mesh.ne,2*mesh.ne);
A1=sparse(2*mesh.ne,2*mesh.ne);
M0=sparse(2*mesh.ne,2*mesh.ne);
M1=sparse(2*mesh.ne,2*mesh.ne);


A0=G0'*D*G0;
A1=G1'*D*G1;
M0=G0_inv*D_inv*G0_inv';
M1=G1_inv*D_inv*G1_inv';


% número de pasos
n=1;
dt=1/(2*n);
M_t0=M0;
M_t1=M1;

for i=1:1
    
% Euler1
M_t1=M_t0-dt*M_t0*(A1-A0)*M_t0;
% Euler2 (atrás)
M_t11=M_t1+dt*M_t1*(A1-A0)*M_t1;
% Promedio Euler
M_t2=(M_t1+M_t11)/2;

% Heun
M_t3=M0-dt/2*(M0*(A1-A0)*M0+M_t0*(A1-A0)*M_t0);
% Heun2 (atrás)
M_t31=M1+dt/2*(M1*(A1-A0)*M1+M_t1*(A1-A0)*M_t1);
% "Promedio Heun"
M_t4=(M_t3+M_t31)/2;

end

%Sin precondicionador
ncond_nopc(j)=cond(full(ABD));

%Nivel 0
ncond_M0(j)=cond(full(M0*ABD));
ncond_prom(j)=cond(full(((M0+M1)/2)*ABD));

%Nivel 1
ncond_euler(j)=cond(full(M_t1*ABD));
ncond_prom_euler(j)=cond(full(M_t2*ABD));

%Nivel 2
ncond_heun(j)=cond(full(M_t3*ABD));
ncond_prom_heun(j)=cond(full(M_t4*ABD));



% Precond. using LU ---- Solución del método Iterativo

% Nivel 0
[xLU0, errorLU0, iterLU0, flagLU0, lambdamaxLU0, condnumberLU0] =ApcgCR2(b*0, b,b, 20000, 0.000001, ABD,M0);
[xLU0p, errorLU0p, iterLU0p, flagLU0p, lambdamaxLU0p, condnumberLU0p] =ApcgCR2(b*0, b,b, 20000, 0.000001, ABD,(M0+M1)/2);

% Nivel 1
[xLU1, errorLU1, iterLU1, flagLU1, lambdamaxLU1, condnumberLU1] =ApcgCR2(b*0, b,b, 20000, 0.000001, ABD,M_t1);
[xLU2, errorLU2, iterLU2, flagLU2, lambdamaxLU2, condnumberLU2] =ApcgCR2(b*0, b,b, 20000, 0.000001, ABD,M_t2);

% Nivel 2
[xLU3, errorLU3, iterLU3, flagLU3, lambdamaxLU3, condnumberLU3] =ApcgCR2(b*0, b,b, 20000, 0.000001, ABD,M_t3);
[xLU4, errorLU4, iterLU4, flagLU4, lambdamaxLU4, condnumberLU4] =ApcgCR2(b*0, b,b, 20000, 0.000001, ABD,M_t4);


Iter0(j)=iterLU0;
Iter0p(j)=iterLU0p;

Iter1(j)=iterLU1;
Iter2(j)=iterLU2;

Iter3(j)=iterLU3;
Iter4(j)=iterLU4;

end


%% Medias número de condición 

% Sin precondicionador
media_ncond_nopc(w)=mean(ncond_nopc);

% Precondicionadores nivel 0
media_ncond_M0(w)=mean(ncond_M0);
media_ncond_prom(w)=mean(ncond_prom);

% Precondicionadores nivel 1
media_ncond_euler(w)=mean(ncond_euler);
media_ncond_prom_euler(w)=mean(ncond_prom_euler);

% Precondicionadores nivel 2
media_ncond_heun(w)=mean(ncond_heun);
media_ncond_prom_heun(w)=mean(ncond_prom_heun);



%% Varianzas número de condición 

% Sin precondicionador
varianza_ncond_nopc(w)=var(ncond_nopc);

% Precondicionadores nivel 0
varianza_ncond_M0(w)=var(ncond_M0);
varianza_ncond_prom(w)=var(ncond_prom);

% Precondicionadores nivel 1
varianza_ncond_euler(w)=var(ncond_euler);
varianza_ncond_prom_euler(w)=var(ncond_prom_euler);

% Precondicionadores nivel 2
varianza_ncond_heun(w)=var(ncond_heun);
varianza_ncond_prom_heun(w)=var(ncond_prom_heun);


%% Medias número de iteraciones 

% Precondicionadores nivel 0
media_iter_M0(w)=mean(Iter0);
media_iter_prom(w)=mean(Iter0p);

% Precondicionadores nivel 1
media_iter_euler(w)=mean(Iter1);
media_iter_prom_euler(w)=mean(Iter2);

% Precondicionadores nivel 2
media_iter_heun(w)=mean(Iter3);
media_iter_prom_heun(w)=mean(Iter4);



%% Varianzas número de iteraciones 

% Precondicionadores nivel 0
varianza_iter_M0(w)=var(Iter0);
varianza_iter_prom(w)=var(Iter0p);

% Precondicionadores nivel 1
varianza_iter_euler(w)=var(Iter1);
varianza_iter_prom_euler(w)=var(Iter2);

% Precondicionadores nivel 2
varianza_iter_heun(w)=var(Iter3);
varianza_iter_prom_heun(w)=var(Iter4);



end

media_ncond_nopc
media_ncond_M0
media_ncond_prom
media_ncond_euler
media_ncond_prom_euler
media_ncond_heun
media_ncond_prom_heun

varianza_ncond_nopc
varianza_ncond_M0
varianza_ncond_prom
varianza_ncond_euler
varianza_ncond_prom_euler
varianza_ncond_heun
varianza_ncond_prom_heun

media_iter_M0
media_iter_prom
media_iter_euler
media_iter_prom_euler
media_iter_heun
media_iter_prom_heun

varianza_iter_M0
varianza_iter_prom
varianza_iter_euler
varianza_iter_prom_euler
varianza_iter_heun
varianza_iter_prom_heun
