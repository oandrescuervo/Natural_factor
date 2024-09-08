function plotvec(vec,MC,M,v,mesh,scale)
ne=mesh.ne;

trimesh(M,v(:,1),v(:,2), v(:,1)*0);
hold on
quiver(MC(:,1),MC(:,2), vec(1:ne),vec((ne+1):2*ne),scale)
view(2)
