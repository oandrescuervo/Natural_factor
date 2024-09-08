clear all
ax=0;bx=1;ay=0;by=1;

n=20; %number of elements in each direction


intx=[ax,bx];inty=[ay,by];

%mesh
[M,B,v,mesh,h,free,MC]=createmeshr(n,intx,inty);
nvel=size(v,1);
[E,DOF,mesh,efree]=edgesenum(M,B,v,mesh,h,free,MC,n);




%coefficient
kappa=coefficient(MC(:,1),MC(:,2));

% Global Neumann matrix for Crouzeixâ€?Raviart (CR) element
r=0.5; % CR mid point based.
[AnCR,bCR]=Nmatrixr(M,DOF,v,mesh,kappa,r);

%Global Neumann matrix of Lagrange (L) elements
[AnL,bL]=NmatrixrP1(M,v,mesh,kappa);


% Diagonal with the cofficient in partial_x and partial_y
D=diag([kappa;kappa]);

%Matrix G_L for partial_x partial_y
[GG1,GG2]=GmatrixP1(M,v,mesh);
list2=1:mesh.nv-1; % we need to take out 1 degree of freedom!
GL=[GG2;-GG1];
clear GG2 GG1
%Matrix G_CR for partial_x partial_y
[G1,G2]=Gmatrixr(M,v,mesh,DOF,r);
GCR=[G1;G2];
clear G1 G2
% epsilon for test

for iii=1:20
    
myeps=10^(10-iii);

%
Gepsilon=[GCR(:,efree),myeps*GL(:,list2)];
Gepsilon2=[GCR(:,efree),(1/myeps)*GL(:,list2)];

G=[GCR(:,efree),GL(:,list2)];

bepsilon=[bCR(efree);myeps*bL(list2)];
b=[bCR(efree);bL(list2)];


% Block digonal matrix 
zaux=sparse(length(list2),mesh.nedgesf);
AL=GL(:,list2)'*D*GL(:,list2);
ACR=GCR(:,efree)'*D*GCR(:,efree);
ABDepsilon=[ACR,zaux';zaux,myeps*AL];
ABD=[ACR,zaux';zaux,myeps*AL];
% Multiplication of 3 square matrices

MatrixM=Gepsilon'*D*Gepsilon2;
%MatrixM=Gepsilon'*D*Gepsilon;
%MatrixM2=Gepsilon'*D*G;
%MatrixM3=G'*D*Gepsilon;

MatrizH1=GCR(:,efree)'*GCR(:,efree);

% free degrees of freedom for block diagonal system
tolfree=[efree,mesh.nedges+(1:length(list2))];

u1=b*0;u2=b*0;u3=b*0;u4=b*0;u5=b*0;


%SOLUCIÓN POR CR INVERTIDO
% CR soluton for comparison
u3=ACR\b(efree);


%SOLUCIÓN CON LA MATRIZ M
% MatrixM solutions using G*D*G
u4(tolfree)=MatrixM\b;

%SOLUCIÓN CON LAS 3 CUADRADAS
% Solution using productof 3 inveres
%u5(tolfree)=G \ ( D \ (G'\b));

% PLOT CR AND MatrixM SOLUTIONS
%     figure
%     subplot(1,2,1)
%     plotscalar(u5,M,DOF,v,mesh)
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

%plot of newumann soltions
%fiugure
% trisurf(M,v(:,1),v(:,2),full((ub))))
% trisurf(M,v(:,1),v(:,2),full((ub2)))



% Preconditionin with M
%PRE=MatrixM;
%PRE=speye(size(ABD));
%[x, error, iter, flag, lambdamax, condnumber] =ApcgCR(b*0, b,b, 20000,0.000001, ABD,PRE);
%u1(tolfree)=x;


% Precond. using LU
[LG,UG]=lu(G);
[xLU, errorLU, iterLU, flagLU, lambdamaxLU, condnumberLU] =ApcgCRLU(b*0, b,b, 20000,0.000001, ABD,LG,UG,D);
u2(tolfree)=xLU;


si=size(u2,1);

i=0;
solCR=zeros(3*n^2-2*n,1);
for kk=1:si
    w=u2(kk);
if abs(w)>0.0001 && i<size(efree,2)
    i=i+1;
    solCR(i)=w;
else
    i;
end
end


solCR2=zeros(3*n^2-2*n,1);
i=0;


    
for kk=1:si
    w=u4(kk);
if abs(w)>0.0001 && i<size(solCR,1)
    i=i+1;
    solCR2(i)=w;
else
    i;
end
end


error(iii,1)=sqrt((u3-solCR2)'*MatrizH1*(u3-solCR2));

end

error