function [G1,G2]=Gmatrixr(M,v,mesh,DOF,r)

nv=mesh.nedges;
ne=mesh.ne;
G1=sparse(ne,nv);
G2=G1;

for i=1:mesh.ne
%   waitbar(i/(2*mesh_f.ne),h)
    %%%%%%%%%%%%%%%%%%  GLOBAL NUMBERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lcol=DOF(i,:); % subdomain indexes of bases.
    %%%%%%%%%%%%%%%%%%% VELOCITY*VLOCITY STIFNESS %%%%%%%%%%%%%%%%%%%%%
    lA=localGr(M(i,:),v,r); % compute local part of A
    G1(i,lcol)=lA(1,:);
    G2(i,lcol)=lA(2,:);
end

