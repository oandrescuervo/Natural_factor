
for i=200:2*mesh.ne
     plotvec(Q(:,i),MC,M,v,mesh,2)
     pause
     hold off
     title(i)
end
