
function [V1,e]=VP(x,nu,l)

M=size(x,1);
C=zeros(M,M);

for i=1:M
    for j=1:M
        % Covarianza Matern
        C(i,j)=covarianza(x(i,:),x(j,:),nu,l);
        % Covarianza cuadrado exponencial 
        %C(i,j)=covarianza_exp(x(i,:),x(j,:));
    end
end

[V,D]=eig(C);
V1=fliplr(V);

e=flipud(eig(C));

