clear all
clc

ax=0;bx=1;ay=0;by=1;

n=20;  %number of elements in each direction
R=100;  %Número de realizaciones MC


intx=[ax,bx];inty=[ay,by];

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
%%%Se puede poner lo que quiera en bL por ser el problema auxiliar



%%%%------------------------------------
%%%% DEPENDECE ON KAPPA
%%%%


% Exponencial con uniforme


% Números aleatorios
rand('seed',1);
va=3*rand(R,4); 
    
for j=1:R
  
% coefficient
kappa=coefficient_contraste(MC(:,1),MC(:,2),va(j,:));




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


for i=1:n
% Euler1
M_t0=M_t0-dt*M_t0*(A1-A0)*M_t0;
% Euler2
M_t1=M_t1+dt*M_t1*(A1-A0)*M_t1;
% "Heun"
M_t2=M0-dt/2*(M0*(A1-A0)*M0+M1*(A1-A0)*M1);
% "Heun 2"
M_t3=M1+dt/2*(M0*(A1-A0)*M0+M1*(A1-A0)*M1);
% "Promedio_euler"
M_t4=(M_t1+M_t0)/2;
% "Promedio_heun"
M_t5=(M_t2+M_t3)/2;
% "Promedio_heun2"
M_t6=(M0+M1)/2;
end

ncond_nopc(j)=cond(ABD);
ncond_euler1(j)=cond(M_t0*ABD);
ncond_euler2(j)=cond(M_t1*ABD);
ncond_heun(j)=cond(M_t2*ABD);
ncond_heun2(j)=cond(M_t3*ABD);
ncond_promE(j)=cond(M_t4*ABD);
ncond_promH(j)=cond(M_t5*ABD);
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
[xLU3, errorLU3, iterLU3, flagLU3, lambdamaxLU3, condnumberLU3] =ApcgCR2(b*0, b,b, 20000, 0.000001, ABD,M_t2);
[xLU4, errorLU4, iterLU4, flagLU4, lambdamaxLU4, condnumberLU4] =ApcgCR2(b*0, b,b, 20000, 0.000001, ABD,M_t5);



%NCond1(j)=condnumberLU;
%NCond2(j)=condnumberLU2;
Iter1(j)=iterLU;
Iter2(j)=iterLU2;
Iter3(j)=iterLU3;
Iter4(j)=iterLU4;

end



media_ncond=mean(ncond_nopc)
media_ncond_euler1=mean(ncond_euler1);
media_ncond_euler2=mean(ncond_euler2)
media_ncond_heun=mean(ncond_heun)
media_ncond_heun2=mean(ncond_heun2)
media_ncond_promE=mean(ncond_promE)
media_ncond_promH=mean(ncond_promH)
media_ncondpc2=mean(ncond_pc2)

var_ncond=var(ncond_nopc)
var_ncond_euler1=var(ncond_euler1)
var_ncond_euler2=var(ncond_euler2)
var_ncond_heun=var(ncond_heun)
var_ncond_heun2=var(ncond_heun2)
var_ncond_promE=var(ncond_promE)
var_ncond_promH=var(ncond_promH)
var_ncondpc2=var(ncond_pc2)


media_iter_promE=mean(Iter2)
media_iter_promH=mean(Iter4)
media_iter_heun=mean(Iter3)
media_iter_pc=mean(Iter1)


var_iter_promE=var(Iter2)
var_iter_pc=var(Iter1)
var_iter_heun=var(Iter3)
var_iter_promH=var(Iter4)