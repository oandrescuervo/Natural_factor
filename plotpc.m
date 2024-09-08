function plotpc(pc,M,v,mesh)

ne=mesh.ne;

for i=1:ne

    trisurf(M(i,:),v(:,1),v(:,2),v(:,1)*0+pc(i))
    hold on
end

view([1,2,1])
