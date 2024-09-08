function plotscalar(p,M,DOF,v,mesh,list)

ne=mesh.ne;

for i=1:ne
    vm=p(DOF(i,:));
    vv=frommidtovertex(M(i,:),vm,v);
    trisurf([1,2,3],v(M(i,:),1),v(M(i,:),2),full(vv),mean(full(vm)))
    hold on
end

view([1,2,1])
