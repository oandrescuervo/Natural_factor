

c=sparse(2*mesh.ne,1)
Iaux=mod(I5,mesh.ne);
c(Iaux)=c(Iaux)+1;
plotpc(c,M,v,mesh)


