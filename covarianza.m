% The Mat�rn class of Covariance functions
% v->infty se obtiene funci�n The squared exponential SE
function cov=covarianza(a,b,v,l)

r=norm(a-b);
if r==0   
    cov=1;   
else 
% gamma: Funci�n gamma 
% besselk(v,k): Modified Bessel function second kind
% v: par�metro, k: entrada de la funci�n
k=(sqrt(2*v)*r/l);
cov=2^(1-v)/gamma(v)*k.^v*besselk(v,k);
%cov=exp(-(norm(a-b))^2/2);
end
