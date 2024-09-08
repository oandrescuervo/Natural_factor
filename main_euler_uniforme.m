clear all
clc

ax=0;bx=1;ay=0;by=1;

n=10;  %number of elements in each direction
R=100;  %Número de realizaciones MC
K=10;  %Número de términos en la serie KL

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

% Números aleatorios
rand('seed',1);
va=1+rand(mesh.ne,R); 
    
for j=1:R
   
% Diagonal with the cofficient in partial_x and partial_y
D=sparse(2*mesh.ne,2*mesh.ne);
D_inv=sparse(2*mesh.ne,2*mesh.ne);
for iaux=1:mesh.ne
    D(iaux,iaux)=va(iaux,j);
    D(mesh.ne+iaux,mesh.ne+iaux)=va(iaux,j);
    D_inv(iaux,iaux)=1/va(iaux,j);
    D_inv(mesh.ne+iaux,mesh.ne+iaux)=1/va(iaux,j);
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

for i=1:n
%Euler1
M_t0=M_t0+dt*M_t0*(A1-A0)*M_t0;
% Euler2
M_t1=M_t1-dt*M_t1*(A1-A0)*M_t1;
% "Heun"
M_t2=M0+dt/2*(M0*(A1-A0)*M0+M1*(A1-A0)*M1);
% "Heun 2"
M_t3=M1-dt/2*(M0*(A1-A0)*M0+M1*(A1-A0)*M1);
% "Promedio"
M_t4=(M0+M1)/2+dt/2*(M0*(A1-A0)*M0-M1*(A1-A0)*M1);
end

ncond_nopc(j)=cond(ABD);
ncond_euler1(j)=cond(M_t0*ABD);
ncond_euler2(j)=cond(M_t1*ABD);
ncond_heun(j)=cond(M_t2*ABD);
ncond_heun2(j)=cond(M_t3*ABD);
ncond_prom(j)=cond(M_t4*ABD);
ncond_pc2(j)=cond(M0*ABD);

% Multiplication of 3 square matrices
%MatrizH1=GCR(:,efree)'*GCR(:,efree);
%MatrixM=G'*D*G;                                %%%Matriz M
%MatrixM2=Gepsilon'*D*G;
%MatrixM3=G'*D*Gepsilon;

% free degrees of freedom for block diagonal system
%tolfree=[efree,mesh.nedges+(1:length(list2))];


%u1=b*0;u2=b*0;u3=b*0;u4=b*0;u5=b*0;


% Precond. using LU ---- Solución del método Iterativo
% Método anterior
[xLU, errorLU, iterLU, flagLU, lambdamaxLU, condnumberLU] =ApcgCR2(b*0, b,b, 20000, 0.000001, ABD,M0);

[xLU2, errorLU2, iterLU2, flagLU2, lambdamaxLU2, condnumberLU2] =ApcgCR2(b*0, b,b, 20000, 0.000001, ABD,M_t4);



%NCond1(j)=condnumberLU;
%NCond2(j)=condnumberLU2;
Iter1(j)=iterLU;
Iter2(j)=iterLU2;




end



media_ncond=mean(ncond_nopc)
media_ncond_euler1=mean(ncond_euler1)
media_ncond_euler2=mean(ncond_euler2)
media_ncond_heun=mean(ncond_heun)
media_ncond_heun2=mean(ncond_heun2)
media_ncond_prom=mean(ncond_prom)
media_ncondpc2=mean(ncond_pc2)

media_iter_prom=mean(Iter2)
media_iter_pc=mean(Iter1)

%media_ncond_prom_aprox=mean(NCond2)
%media_ncondpc2_aprox=mean(NCond1)


