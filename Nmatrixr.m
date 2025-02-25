function [A,b]=Nmatrixr(M,DOF,v,mesh,kappa,r)

% fprintf('Creating Sparse Matrices in subdomain \n');

% Cuadrature points in the reference triangle (See Braess)
%[xi,eta,omega]=setquadrature();


nvel=mesh.nedges; % velocity degrees of freedom==number of vertices.
A=sparse(nvel,nvel);   % grad*grad   
b=sparse(nvel,1);       % right hand side part --->.
%Mul=sparse(npre,1);        % Lagrange multimplier to ensure zero pressure.

% odd loop
% h= waitbar(0,'Please wait...assambling stiffness');
for i=1:mesh.ne
%   waitbar(i/(2*mesh_f.ne),h)
    %%%%%%%%%%%%%%%%%%  GLOBAL NUMBERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lcol=DOF(i,:); % subdomain indexes of bases.
    %%%%%%%%%%%%%%%%%%% VELOCITY*VLOCITY STIFNESS %%%%%%%%%%%%%%%%%%%%%
    lA=localAr(M(i,:),v,r); % compute local part of A
    A(lcol,lcol)=A(lcol,lcol)+kappa(i)*lA;
    %%%%%%%%%%%%%%%%%%%% LOAD VECTOR         %%%%%%%%%%%%%%%%%%%%%%%%%%
    lb=localb(M(i,:),v,mesh,r);
    b(lcol,1)=b(lcol,1)+lb';
    %%%%%%%%%%%%%%%%%%%%%%% LAGRANGE MULTIPLIER %%%%%%%%%%%%%%%%%%%%%%%
    %Mul(M_f(i,:),1)=Mul(M_f(i,:),1)+localmul(M_f(i,:),v_f,mesh_f,xi,eta,omega);
end

