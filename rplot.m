function rplot(vec,MC,M,v,mesh)
ne=mesh.ne;

trimesh(M,v(:,1),v(:,2), v(:,1)*0);
quiver(MC(:,1),MC(:,2), vec(1:ne),MC((ne+1):2*ne))
