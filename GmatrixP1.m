function [G1,G2]=GmatrixP1(M,v,mesh)

nv=mesh.nv;
ne=mesh.ne;
G1=sparse(ne,nv);
G2=G1;

for i=1:mesh.ne
%   waitbar(i/(2*mesh_f.ne),h)
    %%%%%%%%%%%%%%%%%%  GLOBAL NUMBERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lcol=M(i,:); % subdomain indexes of bases.
    %%%%%%%%%%%%%%%%%%% VELOCITY*VLOCITY STIFNESS %%%%%%%%%%%%%%%%%%%%%
    lA=localGP1(M(i,:),v); % compute local part of A
    G1(i,lcol)=lA(1,:);
    G2(i,lcol)=lA(2,:);
end

