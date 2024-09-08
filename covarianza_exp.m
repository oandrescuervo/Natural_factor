% The Matérn class of Covariance functions
% v->infty se obtiene función The squared exponential SE
function cov=covarianza_exp(a,b)
cov=exp(-(norm(a-b))^2/2);
end
