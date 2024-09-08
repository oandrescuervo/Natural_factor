% The Matérn class of Covariance functions
% v->infty se obtiene función The squared exponential SE
function cov=covarianza(a,b,v,l)

r=norm(a-b);
if r==0   
    cov=1;   
else 
% gamma: Función gamma 
% besselk(v,k): Modified Bessel function second kind
% v: parámetro, k: entrada de la función
k=(sqrt(2*v)*r/l);
cov=2^(1-v)/gamma(v)*k.^v*besselk(v,k);
%cov=exp(-(norm(a-b))^2/2);
end
