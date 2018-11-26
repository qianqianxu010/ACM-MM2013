 
load incomp.mat; % the paired comparison data colleted
data=data_ref;
[compare,col]=size(data);
n=16;  % number of nodes
Z=zeros(n,n);
for k = 1:compare
    a=data(k,:);
    Z(a(1),a(2))=Z(a(1),a(2))+1;
end

k=0;
m=sum(sum(Z~=0));
d=zeros(m,n);
w=zeros(1,m);
for i = 1:n
    for j = 1:n
        if (i~=j && Z(i,j)~=0)
            k=k+1;
            w(k)=Z(i,j);
            d(k,i)=1;
            d(k,j)=-1;
        end
    end
end
y=ones(m,1);

index = 1:m;
eps = 1e-10;

[U,S,V] = svd(diag(sqrt(w))*d);
S = diag(S);
r = sum(S>eps);
if r~=n-1
	'Not connected or problematic d'
	return;
end

Sigma = diag(S(1:(n-1)));
V = V(:,1:(n-1));
Uc = U(:,n:m);
Ud = U(:,1:(n-1));
Y = Uc'*diag(sqrt(w))*y;
X = Uc'*diag(sqrt(1./w));

[coefs,lam,nIter] = SolveLasso_FISTA(X,Y);
[ttt,l] = size(coefs);
k_th = l*ones(1,m);
for i = l:-1:1
	k_th(coefs(:,i)~=0)=i;
end
[k_th,order]=sort(k_th);

a=zeros(m,4);
for i = 1:m
	a(i,2)=find(d(order(i),:)==1);
	a(i,3)=find(d(order(i),:)==-1);
end
a(:,1)=lam(k_th);
a(:,4)=w(order);
 
outlier_list=a(:,2:4); % the first and second columns indicate the outlier pair we detected, and the third column is the number of paired comparsion occures on this pair.