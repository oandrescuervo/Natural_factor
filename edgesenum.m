function [E,DOF,mesh,efree]=edgesenum(M,B,v,mesh,h,free,MC,n)
E=sparse(mesh.nv,mesh.nv);
for e=1:mesh.ne
    
    for vi=1:3
        ie=M(e,vi);
        je=M(e,mod(vi,3)+1);
        if ie<=je
            E(ie,je)=1;
        else
            E(je,ie)=1;
        end
    end
    
end
[I,J]=find(E);
EM=[I,J];
mesh.nedges=length(EM);
for ed=1:mesh.nedges
    E(EM(ed,1),EM(ed,2))=ed;
end

E=triu(E)-diag(diag(E))+(triu(E)');

DOF=zeros(size(M));
for e=1:mesh.ne
    d1=E(M(e,1),M(e,2));
    d2=E(M(e,2),M(e,3));
    d3=E(M(e,3),M(e,1));
DOF(e,:)=[d1,d2,d3];
end

eboundary=zeros(4*n,1);
%%
e=0;
%down
for i=1:n
    e=e+1;
    Td=M(B.down(i),:);
    eboundary(e)=E(Td(1),Td(2));
end
%right
for i=1:n
    e=e+1;
    Td=M(B.right(i),:);
    eboundary(e)=E(Td(2),Td(3));
end
%up
for i=1:n
    e=e+1;
    Td=M(B.up(i),:);
    eboundary(e)=E(Td(2),Td(3));
end
%left
for i=1:n
    e=e+1;
    Td=M(B.left(i),:);
    eboundary(e)=E(Td(3),Td(1));
end
efree=setdiff( 1:mesh.nedges, eboundary);
    mesh.nedgesf=length(efree);