
for i=1:mesh.nv
    u(free)=diag(D);
    trisurf(M,v(:,1),v(:,2),full(u))
    pause
    
end
