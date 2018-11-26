function [x,lam,nIter] = SolveLasso_FISTA(A,y,lambda)
if nargin < 3
    lambda_max=norm(A'*y,inf);
    rho=0.01^(1/99);
    lambda=lambda_max*rho.^[0:99];
end
[lam,id]=sort(lambda, 'descend');
[m,n]=size(A);
x=zeros(n,length(lambda));
x0=zeros(n,1);
nIter=[];
warning off
for i=1:length(lam)
    [x_hat,nIter(i), timeSteps, errorSteps] = SolveFISTA(A,y,'initialization',x0,'lambda',lam(i),'stoppingcriterion',3,'tolerance',1e-10);
    x(:,i)=x_hat;
    x0=x_hat;
end
