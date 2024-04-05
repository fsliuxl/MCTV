function Y = diff_y(X,sizeD)
X          = reshape(X,sizeD);
[n1, n2]   = size(X);
Y          = zeros(n1,n2);
Y(:,1)     = X(:,2) - X(:,1);
Y(:,end)   = X(:,end-1) - X(:,end);
for i = 2:n2-1
    Y(:,i) = X(:,i-1) + X(:,i+1) - 2*X(:,i);
end
Y = Y(:);