

for i=1:length(group)
    figure
    vec=diag(D)*0;
 vec(group(i).g)=1;
 plotvec(vec,MC,M,v,mesh,.5)
 axis square
end