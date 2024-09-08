clear all 

ax=0;bx=1;ay=0;by=1;

n=5; %number of elements in each direction
intx=[ax,bx];inty=[ay,by];


for ii=1:4
    n=2*n;
    
    
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
G=[GCR(:,efree),GL(:,list2)];                  %%% Matriz general G

bepsilon=[bCR(efree);myeps*bL(list2)];
b=[bCR(efree);bL(list2)];                      %%% Vector del lado derecho b


%%%Se puede poner lo que quiera en bL por ser el problema auxiliar

%LU factorization of G
[LG,UG]=lu(G);


%%%%------------------------------------
%%%% DEPENDECE ON KAPPA
%%%%


%%%Cambiar esta parte implementando las funciones de covarianza. Tener
%%%cuidado de como generar los vectores propios con los puntos de
%%%cuadratura - ¡¡MIRAR!!

%coefficient
kappa=coefficient(MC(:,1),MC(:,2));
% Diagonal with the cofficient in partial_x and partial_y
D=sparse(2*mesh.ne,2*mesh.ne);
for iaux=1:mesh.ne
    D(iaux,iaux)=kappa(iaux);
    D(mesh.ne+iaux,mesh.ne+iaux)=kappa(iaux);
end

% Block digonal matrix 
zaux=sparse(length(list2),mesh.nedgesf);
AL=GL(:,list2)'*D*GL(:,list2);
ACR=GCR(:,efree)'*D*GCR(:,efree);
ABDepsilon=[ACR,zaux';zaux,myeps*AL];

ABD=[ACR,zaux';zaux,myeps*AL];                  %%%Matriz diagonal por bloques

% Multiplication of 3 square matrices
%MatrixM=G'*D*G;                                %%%Matriz M
%MatrixM2=Gepsilon'*D*G;
%MatrixM3=G'*D*Gepsilon;


MatrizH1=GCR(:,efree)'*GCR(:,efree);

% free degrees of freedom for block diagonal system
tolfree=[efree,mesh.nedges+(1:length(list2))];

u1=b*0;u2=b*0;u3=b*0;u4=b*0;u5=b*0;


% CR soluton for comparison
u3(efree)=ABD(efree,efree)\b(efree);            %%%Solución usando matriz de EF CR



% MatrixM solutions using G*D*G
%u4(tolfree)=MatrixM\b;                           %%%Solución usando la matriz M


% Solution using productof 3 inveres
%u5(tolfree)=G \ ( D \ (G'\b));                   %%% Soluciíon usando la inversa de las cuadradas

%%%Podemos calcular la inversa de D a mano, el computador demora en
%%%reconocer que es diagonal


% PLOT CR AND MatrixM SOLUTIONS
%     figure
%     subplot(1,2,1)
%     plotscalar(u5,M,DOF,v,mesh)                          %%%Gráfica de lasolución de CR (problema original)
%     sss=sprintf("inver3quadradas eps=%.5e",myeps);
%     title(sss)
%     subplot(1,2,2)
%     figure
%     plotscalar(u3,M,DOF,v,mesh)
%     title("CR")
%%



%PLOT DIFFERENCE CR AND MatrixM SOLUTIONS
%diffeps=abs(u3(1:mesh.nedges)-u4(1:mesh.nedges));
%figure
%plotscalar(diffeps,M,DOF,v,mesh)
%title("difference")
%colorbar
%view(2)
%axis square

% SOLUTIONS OF NEUMANN PROBLEM COMING FROM 
% THE BLOCK SYSTEM MatrixM
ub=bL*0;
%ub(list2)=u5(mesh.nedges+(1:length(list2)));
%ub=ub-mean(ub);

% P1 SOLUTION OF NEUMANN PROBLEM
ub2=bL*0;
%AnLf=AnL(list2,list2);
%ub2(list2)=AnLf\bL(list2);
ub2=ub2-mean(ub2);

%plot of Neumann solutions
%figure
% trisurf(M,v(:,1),v(:,2),full((ub))))
%trisurf(M,v(:,1),v(:,2),full((ub2)))                  %%%Gráfica de la solución con Neumann (auxiliar)




% Preconditionin with M
%PRE=MatrixM;
%PRE=speye(size(ABD));
%[x, error, iter, flag, lambdamax, condnumber] =ApcgCR(b*0, b,b, 20000,0.000001, ABD,PRE);
%u1(tolfree)=x;


% Precond. using LU ---- Solución del método Iterativo
%%%No hay necesidad de montar toda la matriz (mirar al interior de la función), solo GL, GR y D.
[xLU, errorLU, iterLU, flagLU, lambdamaxLU, condnumberLU] =ApcgCRLU(b*0, b,b, 2000000,0.000001, ABD,LG,UG,D);
u2(tolfree)=xLU;

[valx,valy]=midpoint(v,E,efree);

si=size(u2,1);


exact=analytic(valx,valy);


i=0;
solCR=zeros(3*n^2-2*n,1);
for kk=1:si
    w=u2(kk);
if abs(w)>0.000001
    i=i+1;
    solCR(i)=w;
else
    i;
end
end

Error(ii,1)=sqrt((exact-solCR)'*MatrizH1*(exact-solCR));

end
Error
